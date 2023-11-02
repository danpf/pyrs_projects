from typing import cast, Iterator, Any
import requests
import tarfile
import io
import csv
from tqdm import tqdm
from pathlib import Path
import asyncio
from time import sleep

import aiohttp
from sqlalchemy import create_engine, not_, func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.sqlite import insert

from .models import CombinedData, TraceData, Base, UrlQueryData, DataPackageData, TraceRawData


DATABASE_PATH = "sqlite:///boldsystems.db"

engine = create_engine(DATABASE_PATH)
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)


def parse_tsv_weird_format(input_str: str) -> list[dict[str, Any]]:
    lines = input_str.split("\n")
    headers = lines[0].split("\t")

    data = []
    current_row = []

    for line in lines[1:]:
        fields = line.split("\t")

        if len(fields) == len(headers):
            data.append(dict(zip(headers, fields)))
        elif len(fields) > len(headers):
            raise RuntimeError(f"Too many fields? {len(fields)=} {len(headers)=}")
        else:
            if not current_row and len(fields) < len(headers):
                current_row = fields
            else:
                current_row[-1] += "\n" + fields[0]
                current_row += fields[1:]
                if len(current_row) == len(headers):
                    data.append(dict(zip(headers, current_row)))
                    current_row = []
                elif len(current_row) > len(headers):
                    raise RuntimeError(f"Too many fields? {len(fields)=} {len(headers)=}")
    return data


def tsv_to_dict_iterator(tsv_file: Path) -> Iterator[dict[str, Any]]:
    with open(tsv_file, "r", encoding="utf-8", errors="ignore") as f:
        header = next(f)
        header_d = header.strip().split("\t")
        for line in f:
            if line:
                cline = line.strip().split("\t")
                crow = dict(zip(header_d, cline))
                yield crow


async def download_data(url: str, session: Any | None = None, timeout: int | None = None) -> bytes:
    if session is None:
        session = Session()
    cached_data = session.query(UrlQueryData).filter_by(url=url).first()
    if cached_data is not None:
        print("found cache?")
        return cast(bytes, cached_data.url_response)
    async with aiohttp.ClientSession() as aio_session:
        async with aio_session.get(url, timeout=timeout) as response:
            data_o = io.BytesIO()
            async for data in response.content.iter_chunked(1024):
                data_o.write(data)
            if response.status == 200:
                datavalue = data_o.getvalue()
                cache_entry = UrlQueryData(url=url, url_response=datavalue)
                session.add(cache_entry)
                session.commit()
                return datavalue
            else:
                raise RuntimeError(f"Error: Unable to download data. Status code: {response.status} {url}")


async def download_and_process_bytes(sesh, url):
    tries = 3
    while True:
        try:
            tsv_bytes = await download_data(url, sesh, 60)
            break
        except:
            tries -= 1
            if tries == 0:
                print("Problem with url?", url)
                raise
    rows = parse_tsv_weird_format(io.StringIO(tsv_bytes.decode(encoding="utf-8", errors="ignore")).getvalue())

    for row in rows:
        if not row["trace_ids"]:
            row["trace_ids"] = "_"
        stmt = (
            insert(CombinedData)
            .values(**row)
            .on_conflict_do_nothing(index_elements=["processid", "sampleid", "recordID", "trace_ids"])
        )
        sesh.execute(stmt)
    sesh.commit()


async def download_and_process_trace_data(session, trace_id, trace_url, trace_name, direction, seq_primer, marker_code, processid):
    cached_data = session.query(TraceRawData).filter_by(trace_id=trace_id, trace_url_id=trace_url).first()
    if cached_data is not None:
        return cached_data.tracebytes
    while True:
        try:
            async with aiohttp.ClientSession() as aio_session:
                print("downloading", trace_url, processid)
                async with aio_session.get(trace_url, timeout=20) as response:
                    data_o = io.BytesIO()
                    async for data in response.content.iter_chunked(1024):
                        data_o.write(data)
                    if response.status == 200:
                        datavalue = data_o.getvalue()
                        cache_entry = TraceRawData(
                            trace_id=trace_id,
                            trace_url_id=trace_url,
                            trace_name=trace_name,
                            tracebytes=datavalue,
                            direction=direction,
                            seq_primer=seq_primer,
                            marker_code=marker_code,
                            processid=processid,
                        )
                        session.add(cache_entry)
                        session.commit()
                        return datavalue
                    else:
                        raise RuntimeError(f"Error: Unable to download data. Status code: {response.status} {trace_url}")
        except Exception as e:
            print(f"Caught exception {e} from {trace_url=} so sleeping...")
            await asyncio.sleep(60)



class BoldSystemsAPI:
    def __init__(self):
        self.base_combined_url = "http://www.boldsystems.org/index.php/API_Public/combined?"
        self.base_trace_url = "http://www.boldsystems.org/index.php/API_Public/trace?"

    def build_query(self, base_url, params):
        query = base_url
        to_add = []
        for key, value in params.items():
            if value:
                to_add.append(f"{key}={value}")
        query += "&".join(to_add)
        return query

    def get_data_package_data(self, data_package_file: Path):
        session = Session()
        count = 10_000
        for row in tsv_to_dict_iterator(data_package_file):
            row["class_"] = row.pop("class")
            try:
                stmt = insert(DataPackageData).values(**row).on_conflict_do_nothing(index_elements=["processid"])
            except:
                print("fail", row)
                continue
            session.execute(stmt)
            count -= 1
            if count == 0:
                session.commit()
                count = 10_000
        session.commit()

    async def get_combined_data_from_local_package_data(self):
        session = Session()

        params = {
            "format": "tsv",
        }
        url_no_ids = self.build_query(self.base_combined_url, params)
        base_url_len = len(url_no_ids)

        def get_combined_str(p_ids: list[str]) -> str:
            return "|".join(p_ids)

        c_pids = []
        tasks = []
        for row in tqdm(
            session.query(DataPackageData)
            .filter(not_(DataPackageData.processid.in_(session.query(CombinedData.processid))))
            .yield_per(1000),
            total=9_100_000,
        ):
            c_pids.append(str(row.processid))
            combo_str = get_combined_str(c_pids)
            if len(combo_str) + base_url_len > 2000:
                params.update({"ids": combo_str})
                url = self.build_query(self.base_combined_url, params)
                # download_and_process_bytes(session, url)
                task = asyncio.create_task(download_and_process_bytes(session, url))
                tasks.append(task)
                if len(tasks) >= 20:
                    await asyncio.gather(*tasks)
                    tasks = []
                c_pids = []
        if c_pids:
            params.update({"ids": get_combined_str(c_pids)})
            url = self.build_query(self.base_combined_url, params)
            task = asyncio.create_task(download_and_process_bytes(session, url))
            tasks.append(task)
            if len(tasks) >= 20:
                await asyncio.gather(*tasks)

    def get_combined_data(self, search_params: dict[str, str]):
        params = {
            "taxon": search_params.get("taxon"),
            "ids": search_params.get("ids"),
            "bin": search_params.get("bin"),
            "container": search_params.get("container"),
            "institutions": search_params.get("institutions"),
            "researchers": search_params.get("researchers"),
            "geo": search_params.get("geo"),
            "marker": search_params.get("marker"),
            "format": "tsv",
        }

        url = self.build_query(self.base_combined_url, params)
        tsv_bytes = download_data(url)
        reader = csv.reader(io.StringIO(tsv_bytes.decode(encoding="utf-8", errors="ignore")), delimiter="\t")

        session = Session()
        try:
            header = next(reader)
        except StopIteration:
            print(f"Problem with call: {url=} {tsv_bytes=}")

        # seen = set()
        for row in reader:
            ddd = dict(zip(header, row))
            # dd = ddd["processid"]
            # print(ddd)
            # if dd in seen:
            #     if ddd["trace_ids"]:
            #         print("skipping dupe", dd, ddd)
            #     continue
            # else:
            #     seen.add(dd)
            combined_data = CombinedData(**ddd)
            if not cast(str, combined_data.trace_ids):
                continue
            if (
                session.query(CombinedData)
                .filter_by(
                    processid=combined_data.processid,
                    sampleid=combined_data.sampleid,
                    recordID=combined_data.recordID,
                    trace_ids=combined_data.trace_ids,
                )
                .first()
                is None
            ):
                session.add(combined_data)
                session.commit()

    async def download_all_trace_data(self, session):
        tasks = []
        for row in session.query(CombinedData).filter(func.length(CombinedData.trace_ids) > 1).yield_per(1000):
            traces = row.trace_ids.split("|")
            trace_urls = row.trace_links.split("|")
            trace_names = row.trace_names.split("|")
            directions = row.directions.split("|")
            seq_primers = row.seq_primers.split("|")
            marker_codes = row.marker_codes.split("|")

            for trace, trace_url, trace_name, direction, seq_primer, marker_code in zip(
                traces, trace_urls, trace_names, directions, seq_primers, marker_codes, strict=True
            ):
                if trace_url.endswith("/"):
                    # sometimes ids and urls are there, but the url doesn't have a download link, it just ends after the /
                    continue
                tasks.append(
                    asyncio.create_task(
                        download_and_process_trace_data(session, trace, trace_url, trace_name, direction, seq_primer, marker_code, row.processid)
                    )
                )
                if len(tasks) == 5:
                    await asyncio.gather(*tasks)
                    tasks = []
        if tasks:
            await asyncio.gather(*tasks)
            tasks = []

        session.commit()

    def get_trace_data(self, search_params, session):
        params = {
            "taxon": search_params.get("taxon"),
            "ids": search_params.get("ids"),
            "bin": search_params.get("bin"),
            "container": search_params.get("container"),
            "institutions": search_params.get("institutions"),
            "researchers": search_params.get("researchers"),
            "geo": search_params.get("geo"),
            "marker": search_params.get("marker"),
        }

        query_url = self.build_query(self.base_trace_url, params)
        tar_bytes = download_data(query_url)
        # tar_bytes = open("mammalia_canada_trace.tar", "rb").read()
        with tarfile.open(name=None, fileobj=io.BytesIO(tar_bytes)) as tar:
            trace_info_file = tar.getmember("TRACE_FILE_INFO.txt")
            trace_info_file_io = tar.extractfile(trace_info_file)
            if trace_info_file_io is None:
                raise RuntimeError(f"No result from {query_url}")
            else:
                reader = csv.reader(io.StringIO(trace_info_file_io.read().decode()), delimiter="\t")

            header = next(reader)  # Skip the header row

            count = 0
            for row in reader:
                d_row = dict(zip(header, row))
                if not d_row:
                    continue
                rawfilename = d_row["TRACEFILE"].split("/")[-1].split("+")[0] + ".ab1"
                tracefile_member = tar.getmember(rawfilename)
                tracebytes_io = tar.extractfile(tracefile_member)
                if tracebytes_io is None:
                    continue
                else:
                    tracebytes = tracebytes_io.read()

                # Create an instance of TraceData and populate its attributes
                ov: dict[str, str | bytes] = {
                    "PROCESSID": "CCSMA054-07",
                    "TAXON": "Hemithiris psittacea",
                    "MARKER": "COI-5P",
                    "GENBANK_ACCESSION": "",
                    "TRACEFILE": "CCSMA/CCSMA054-07[LCO1490_t1,HCO2198_t1]_R+1248800947.ab1",
                }
                ov = {k.lower(): v for k, v in d_row.items()}
                ov["rawfilename"] = rawfilename
                ov["tracebytes"] = tracebytes
                entry = TraceData(**ov)
                # Add the entry to the session and commit
                session.add(entry)
                count += 1
            if count:
                session.commit()
