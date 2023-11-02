

from dpf_sanger_sequencing.boldsystems.bold_api_access import BoldSystemsAPI

api = BoldSystemsAPI()
search_params = {"taxon": "Brachiopoda"}
# api.get_combined_data(search_params, "tsv")
api.get_trace_data(search_params)
