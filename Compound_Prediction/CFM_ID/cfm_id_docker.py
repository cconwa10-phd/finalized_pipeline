import os
import subprocess

os.chdir("/Users/ciaraconway")
#print(os.getcwd())
#smile = "O=C(O)CCC1=NC(C=2C=CC=CC2)=C(O1)C=3C=CC=CC3"
#print(smile)
#result = os.system("docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cfm-predict " + "'" + smile + "'" + " 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1\"")
#print(result)



def cfm_id_predict_loop(list):
    for smile in list:
        os.system("docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cfm-predict " + "'" + smile + "'" + " 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1\"")




cmpd_list = ["O=C(O)CCC1=NC(C=2C=CC=CC2)=C(O1)C=3C=CC=CC3", "CC1=NC=C(N1CC(CCl)O)[N+](=O)[O-]", "N=C1N=C(O)C(O1)C=2C=CC=CC2", "O=C1N=C(O)C(=O)N1", "O=P(O)(O)OCCN"]

cfm_id_predict_loop(cmpd_list)
