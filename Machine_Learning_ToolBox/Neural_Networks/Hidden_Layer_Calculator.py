import pandas as pd

list_list = list()
alphas = [2, 6, 8, 10]
features = [5, 10, 15, 20, 30]
for alpha in alphas:
    for feature in features:
        num_sam = 569
        num_inputs = feature
        num_outputs = 2
        hidd = round(num_sam/(alpha*(num_inputs+num_outputs)))
        list_list.append((alpha, feature, hidd))
df = pd.DataFrame(list_list, columns=["alpha", "feature", "# neurons in hidden layer"])
df.to_csv("/Users/ciaraconway/Documents/numhiddlayer.csv")

