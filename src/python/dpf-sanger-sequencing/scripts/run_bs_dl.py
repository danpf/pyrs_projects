from pathlib import Path
from dpf_sanger_sequencing.boldsystems import BoldSystemsAPI, Session, CombinedData
from dpf_sanger_sequencing.boldsystems.bold_api_access import parse_tsv_weird_format
import asyncio


def store_combined_data_for_taxonomies(taxonomies: list[str]):
    api = BoldSystemsAPI()
    # api.get_data_package_data(Path("~/Downloads/bold_pub/BOLD_Public.05-May-2023.tsv").expanduser())
    # print("x1?")
    # asyncio.run(api.get_combined_data_from_local_package_data())
    # raise

    # for i, taxonomy in enumerate(taxonomies):
    #     print(f"getting {i} / {len(taxonomies)} {taxonomy=}")
    #     api.get_combined_data(dict(taxon=taxonomy))
    session = Session()
    # rows = session.query(CombinedData).all()
    asyncio.run(api.download_all_trace_data(session))

    # for row in rows:
        # print(row.processid)


if __name__ == "__main__":
    taxonomies = [
        "Acanthocephala",
        # "Acoelomorpha",  # Missing data?
        "Annelida",
        "Arthropoda",
        "Brachiopoda",
        "Bryozoa",
        "Chaetognatha",
        "Chordata",
        "Cnidaria",
        "Ctenophora",
        "Cycliophora",
        "Echinodermata",
        "Entoprocta",
        "Gastrotricha",
        "Gnathostomulida",
        "Hemichordata",
        "Kinorhyncha",
        "Mollusca",
        "Nematoda",
        "Nematomorpha",
        "Nemertea",
        "Onychophora",
        "Phoronida",
        "Placozoa",
        "Platyhelminthes",
        "Porifera",
        "Priapulida",
        "Rhombozoa",
        "Rotifera",
        "Tardigrada",
        "Xenacoelomorpha",
        "Bryophyta",
        "Chlorophyta",
        "Lycopodiophyta",
        "Magnoliophyta",
        "Pinophyta",
        "Pteridophyta",
        "Ascomycota",
        "Basidiomycota",
        "Chytridiomycota",
        "Glomeromycota",
        "Myxomycota",
        "Zygomycota",
        "Chlorarachniophyta",
        "Ciliophora",
        "Heterokontophyta",
        "Pyrrophycophyta",
        "Rhodophyta",
    ]
    store_combined_data_for_taxonomies(taxonomies)
