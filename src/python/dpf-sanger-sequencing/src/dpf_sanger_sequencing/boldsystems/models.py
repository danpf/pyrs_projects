from sqlalchemy import Column, String, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, ForeignKey, create_engine, Table
from sqlalchemy.orm import relationship, backref, Mapped, mapped_column
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class CombinedData(Base):
    __tablename__ = "combined_data"

    processid = Column(String, primary_key=True)
    sampleid = Column(String, primary_key=True)
    recordID = Column(String, primary_key=True)
    catalognum = Column(String)
    fieldnum = Column(String)
    institution_storing = Column(String)
    collection_code = Column(String)
    bin_uri = Column(String)
    phylum_taxID = Column(String)
    phylum_name = Column(String)
    class_taxID = Column(String)
    class_name = Column(String)
    order_taxID = Column(String)
    order_name = Column(String)
    family_taxID = Column(String)
    family_name = Column(String)
    subfamily_taxID = Column(String)
    subfamily_name = Column(String)
    genus_taxID = Column(String)
    genus_name = Column(String)
    species_taxID = Column(String)
    species_name = Column(String)
    subspecies_taxID = Column(String)
    subspecies_name = Column(String)
    identification_provided_by = Column(String)
    identification_method = Column(String)
    identification_reference = Column(String)
    tax_note = Column(String)
    voucher_status = Column(String)
    tissue_type = Column(String)
    collection_event_id = Column(String)
    collectors = Column(String)
    collectiondate_start = Column(String)
    collectiondate_end = Column(String)
    collectiontime = Column(String)
    collection_note = Column(String)
    site_code = Column(String)
    sampling_protocol = Column(String)
    lifestage = Column(String)
    sex = Column(String)
    reproduction = Column(String)
    habitat = Column(String)
    associated_specimens = Column(String)
    associated_taxa = Column(String)
    extrainfo = Column(String)
    notes = Column(String)
    lat = Column(String)
    lon = Column(String)
    coord_source = Column(String)
    coord_accuracy = Column(String)
    elev = Column(String)
    depth = Column(String)
    elev_accuracy = Column(String)
    depth_accuracy = Column(String)
    country = Column(String)
    province_state = Column(String)
    region = Column(String)
    sector = Column(String)
    exactsite = Column(String)
    image_ids = Column(String)
    image_urls = Column(String)
    media_descriptors = Column(String)
    captions = Column(String)
    copyright_holders = Column(String)
    copyright_years = Column(String)
    copyright_licenses = Column(String)
    copyright_institutions = Column(String)
    photographers = Column(String)
    sequenceID = Column(String)
    markercode = Column(String)
    genbank_accession = Column(String)
    nucleotides = Column(String)
    trace_ids = Column(String, primary_key=True)
    trace_names = Column(String)
    trace_links = Column(String)
    run_dates = Column(String)
    sequencing_centers = Column(String)
    directions = Column(String)
    seq_primers = Column(String)
    marker_codes = Column(String)


class DataPackageData(Base):
    __tablename__ = "data_package_data"

    processid = Column(String, primary_key=True)
    sampleid = Column(String)
    specimenid = Column(String)
    museumid = Column(String)
    fieldid = Column(String)
    inst = Column(String)
    bin_uri = Column(String)
    identification = Column(String)
    funding_src = Column(String)
    kingdom = Column(String)
    phylum = Column(String)
    class_ = Column(String)
    order = Column(String)
    family = Column(String)
    subfamily = Column(String)
    genus = Column(String)
    species = Column(String)
    subspecies = Column(String)
    identified_by = Column(String)
    voucher_type = Column(String)
    collectors = Column(String)
    collection_date = Column(String)
    collection_date_accuracy = Column(String)
    life_stage = Column(String)
    sex = Column(String)
    reproduction = Column(String)
    extrainfo = Column(String)
    notes = Column(String)
    coord = Column(String)
    coord_source = Column(String)
    coord_accuracy = Column(String)
    elev = Column(String)
    depth = Column(String)
    elev_accuracy = Column(String)
    depth_accuracy = Column(String)
    country = Column(String)
    province = Column(String)
    country_iso = Column(String)
    region = Column(String)
    sector = Column(String)
    site = Column(String)
    collection_time = Column(String)
    habitat = Column(String)
    collection_note = Column(String)
    associated_taxa = Column(String)
    associated_specimen = Column(String)
    species_reference = Column(String)
    identification_method = Column(String)
    recordset_code_arr = Column(String)
    gb_acs = Column(String)
    marker_code = Column(String)
    nucraw = Column(String)
    sequence_run_site = Column(String)
    processid_minted_date = Column(String)
    sequence_upload_date = Column(String)
    identification_rank = Column(String)


class TraceData(Base):
    __tablename__ = "trace_data"

    processid = Column(String)
    taxon = Column(String)
    marker = Column(String)
    genbank_accession = Column(String)
    tracefile = Column(String, primary_key=True)
    rawfilename = Column(String)
    tracebytes = Column(LargeBinary)


class UrlQueryData(Base):
    __tablename__ = "urlquery_data"

    url = Column(String, primary_key=True)
    url_response = Column(LargeBinary)


class TraceRawData(Base):
    __tablename__ = "raw_traces"

    trace_id = Column(String, primary_key=True)
    trace_url_id = Column(String, primary_key=True)
    processid = Column(String, primary_key=True)
    tracebytes = Column(LargeBinary)
    direction = Column(String)
    trace_name = Column(String)
    seq_primer = Column(String)
    marker_code = Column(String)

    @staticmethod
    def base_url_id_url():
        return "http://trace.boldsystems.org/traceIO/bold.org/"


class TrainingData(Base):
    __tablename__ = "training_data"

    trace_names: Mapped[str] = mapped_column(primary_key=True)
    target_sequence: Mapped[str] = ForeignKey("table1.id")
    mapped_column()


TrainingData = Table(
    "TrainingData",
    Base.metadata,
    Column("table1_id", Integer, ForeignKey("table1.id")),
    Column("table2_id", Integer, ForeignKey("table2.id")),
)

Table1.associated_rows = relationship(
    "Table2", secondary=association_table, backref=backref("table1_rows", lazy="dynamic"), lazy="dynamic"
)

Table2.associated_rows = relationship(
    "Table1", secondary=association_table, backref=backref("table2_rows", lazy="dynamic"), lazy="dynamic"
)
