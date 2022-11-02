#creator: Ciara Conway

import pandas as pd
import numpy as np
import json
import requests
import collections


### Matches HMDB Ids with KEGG Ids using pandas merge - Need to have hmdb csv and phytochemical spreadsheet
def hmdb_xml_match_merge(target_df, ref_df, file):
    target_left = pd.read_csv(target_df, dtype=str)
    ref_right = pd.read_csv(ref_df, dtype=str)
    print(target_left.head(), ref_right.head())
    hmdb_matches = pd.merge(target_left, ref_right, how="left", on=["HMDB"])
    hmdb_matches = hmdb_matches.fillna("Nah")
    hmdb_matches.to_csv(file)
    print(hmdb_matches.head())
    return hmdb_matches

# Calls KEGG with matching Kegg IDs to find pathway Ids
def kegg_pull(cmpd_kegg_id_df):
    map_dict = collections.defaultdict(list)
    kegg_df = pd.read_csv(cmpd_kegg_id_df, dtype=str)
    ids = kegg_df["kegg_id"]
    for id in ids:
        print(id)
        if id.startswith("C"):
            urlpattern_kegg = 'http://rest.kegg.jp/link/pathway/' + id
            response = requests.get(urlpattern_kegg)
            data = response.text
            data = data.strip().split('\n')
            for i in data:
                if i == '':
                    place_holder = id
                    map_dict[place_holder] = ["None"]
                else:
                    k, v = i.split('\t', 1)
                    k = k.split(':', 1)[1]
                    v = v.split(':', 1)[1]
                    #kegg_path(id) ### Call image acquisition here
                    if k in map_dict.keys():
                        map_dict[k].append(v)
                    else:
                        map_dict[k] = [v]
        else:
            pass
    df = {"ID":ids, "KEGG_ID_Url":map_dict.keys(), "Pathways":map_dict.values()}
    dataframe = pd.DataFrame(df, columns= ["ID", "KEGG_ID_Url", "Pathways"])
    dataframe.to_csv("/Users/ciaraconway/Documents/Kay_Lab/kegg_test_result.csv")
    return dataframe

# Calls KEGG to Generate image of pathway
def kegg_path(id):
    urlpattern_map = 'http://rest.kegg.jp/get/' + id + '/image'
    response = requests.get(urlpattern_map)
    file = open('test_image' + id + '.png', 'wb')
    file.write(response.content)
    file.close()

### Calls reactome with matching Kegg IDs to find Dbids, Stids, pathway Ids
def reactome_pull(hmdb_matches_df):
    pathway_dict = dict()
    ids = hmdb_matches_df["kegg_id"]
    dbid_list = []
    db_pathway_list = []
    pathway_list_df = []
    for value in ids:
        if value.startswith("C"): #Only loops through those values that have a KEGG ID
            db_id_string = ""
            pathway_string = ""
            bare_string = ""
            id_reactome_list = []
            url_reactome_xreg = 'https://reactome.org/ContentService/references/mapping/' + value
            headers = {
            'accept': 'application/json',
            }
            response = requests.get(url_reactome_xreg,headers=headers) #Requests Dbid
            json_file_id = response.text
            code_response = response.status_code
            id_return = json.loads(json_file_id)
            count = 0
            for i in id_return:
                if code_response == 404:
                    while count == 0:
                        db_id_string += value
                        pathway_string += value
                        bare_string += "Nah"
                        count += 1
                else:
                    print(i["dbId"])
                    dbid = i["dbId"]
                    id_reactome_list.append(i["dbId"])
                    db_id_string += str(i["dbId"]) + "\n"
                    url_reactome_id_stid = 'https://reactome.org/ContentService/data/query/enhanced/' + str(i["dbId"])
                    response_enhanced = requests.get(url_reactome_id_stid,headers=headers) #Requests Stid through enhanced search
                    code_response_enhanced = response_enhanced.status_code
                    json_file_enhanced = response_enhanced.text
                    stid_return = json.loads(json_file_enhanced)
                    stid_list = []
                    count1 = 0
                    for v in stid_return['physicalEntity']: #Loops through all Stids, multiple Stids map to single Dbid
                        if code_response_enhanced == 404:
                            while count1 == 0:
                                pathway_string += "Nah"
                                bare_string += "Nah"
                                count1 += 1
                        else:
                            pathway_list = []
                            stid = v['stId']
                            stid_list.append(stid)
                            pathway_string += "Parent_stID: " + stid + "\n"
                            bare_string += "Parent_stID: " + stid + "\n"
                            url_reactome_stid_path = 'https://reactome.org/ContentService/data/pathways/low/entity/' + str(stid)
                            response_pathway = requests.get(url_reactome_stid_path,headers=headers) #Requests pathway information from Stid
                            code_response_pathway = response_pathway.status_code
                            json_file_path = response_pathway.text
                            pathway_return = json.loads(json_file_path)
                            count2 = 0
                            for pathway in pathway_return:
                                if code_response_pathway == 404:
                                    while count2 == 0:
                                        pathway_string += "Nah"
                                        bare_string += "Nah"
                                        count2 += 1
                                elif pathway['speciesName'] == 'Homo sapiens':
                                    pathway_list.append((pathway["displayName"], pathway["dbId"], pathway["stId"], pathway["isInDisease"]))
                                    pathway_string += "Name: " + pathway["displayName"] + "\n" + "dbId_pathway: " + str(pathway["dbId"]) + "\n" + "stId_pathway: " + pathway["stId"] + "\n" + "Associated_w/_Disease?: " + str(pathway["isInDisease"]) + "\n"
                                    bare_string += pathway["displayName"] + ":" + pathway["stId"] + "\n"
                                    pathway_list_df.append((value, dbid, stid, pathway["displayName"], pathway["stId"], pathway["isInDisease"]))
                                    pathway_dict[pathway["stId"]] = pathway["displayName"]

            print(pathway_string)
            print(bare_string)
            db_pathway_list.append((db_id_string, pathway_string, bare_string))
            dbid_list.append(db_id_string)
        else:
            db_pathway_list.append(("Nah", "Nah", "Nah"))
            dbid_list.append("Nah")
    print(db_pathway_list)
    #Separate dataframe for pathway information only
    pathwaydf = pd.DataFrame(pathway_list_df, columns=["value", "db_id", "stid","pathway", "pathwayid", "associated with disease"])
    pathwaydf.to_csv("/Users/cconwa10/Desktop/reactome_ids_pathonly.csv")
    id_df = pd.DataFrame(db_pathway_list, columns=["db_id", "pathway_info", "pathway_bare"])
    #Separate dataframe for pathway info linked to Dbid
    id_df.to_csv("/Users/cconwa10/Desktop/reactome_ids.csv")
    hmdb_matches_df["dbid"] = dbid_list
    #Separate dataframe for hmdb merged file with associated Dbids
    hmdb_matches_df.to_csv("/Users/cconwa10/Desktop/hmdbmerge.csv")
    hmdb_matches_df = pd.concat([hmdb_matches_df, id_df], axis=1)
    #Separate dataframe for concatenated hmdb df and all associated reactome info
    #Note: This dataframe is not reliable for all entries since some require more characters than excel cells will allow
    hmdb_matches_df.to_csv("/Users/cconwa10/Desktop/reactome_ids_full.csv")
    print(hmdb_matches_df.head(10))
    return hmdb_matches_df


# Gather Pathway JPEGs Reactome
def reactomeJpeg(file):
    headers = {
        'accept': 'image/png',
    }

    pathwayDF = pd.read_csv(file)
    pathwayIDs = pathwayDF["pathwayid"].values

    for i, pathwayID in enumerate(pathwayIDs):
        flagId = pathwayDF.loc[pathwayDF.pathwayid == pathwayIDs, "stid"][i]
        name = pathwayDF.loc[pathwayDF.pathwayid == pathwayIDs, "pathway"][i]
        keggId = pathwayDF.loc[pathwayDF.pathwayid == pathwayIDs, "value"][i]

        urlJPEG = "https://reactome.org/ContentService/exporter/diagram/" + pathwayID +".jpeg?quality=10&flg=" + flagId + "&flgInteractors=true&title=true&margin=15&ehld=true&diagramProfile=Modern&resource=TOTAL&analysisProfile=Standard"
        response = requests.get(urlJPEG, headers=headers, allow_redirects=True)
        codeResponse = response.status_code

        if codeResponse == 200:
            path = "/Users/cconwa10/Documents/PR_JPEG/"
            fileName = keggId + ":" + flagId + ":" + pathwayID + ".png"
            pngFile = open(path + fileName, "wb")
            pngFile.write(response.content)
            pngFile.close()
            print(pathwayID)
        else:
            pass




def main():
    pass


if __name__ == '__main__':
    main()
