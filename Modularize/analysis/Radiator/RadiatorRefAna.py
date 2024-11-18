"""
This program focus on analyze the references.
"""
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', "..", ".."))
from Modularize.analysis.Radiator.RadiatorSetAna import main_analysis, get_references_from_ResultJson, save_ref


""" Manually fill in """
target_q = 'q0'
conditional_folder = "Radiator_WS"          # the previous folder from temperature folder     
sample_folder = "Radiator_wisconsinQ1"     # the previous folder from conditional_folder

# ? Save the references from a mK set
ref_info = {"before":{"save_ref_to_json":False, "mK_folder_path":""},"recover":{"save_ref_to_json":False, "mK_folder_path":""}}

# save reference
for ref_type in ref_info:
    if ref_info[ref_type]["save_ref_to_json"]:
        if ref_info[ref_type]["mK_folder_path"] != "":
            main_analysis(target_q, ref_info[ref_type]["mK_folder_path"])
            ref_dict = get_references_from_ResultJson(ref_info[ref_type]["mK_folder_path"])
            save_ref(ref_dict,target_q,sample_name=sample_folder,conditional_name=conditional_folder,ref_type=ref_type)