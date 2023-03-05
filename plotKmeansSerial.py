import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

arg=sys.argv

if len (arg) < 3 : 
    print("Missing arguments : <csv_file_name> <log_directory> <no_of_clustures> <output_file_name>")
    exit()

else :
    csv_file = arg[1]
    log_dir = arg[2]
    nclusters = int(arg[3])
    output_file_name = arg[4]

plt.figure()

df = pd.read_csv(log_dir + '/' + csv_file)
sns.scatterplot(x=df.x, y=df.y, 
                hue=df.cluster_id, 
                palette=sns.color_palette("hls", n_colors=nclusters))
plt.xlabel("Annual income (k$)")
plt.ylabel("Spending Score (1-100)")
plt.title("Clustered: spending (y) vs income (x)")

plt.savefig(log_dir +'/'+output_file_name +'_serial.png')

