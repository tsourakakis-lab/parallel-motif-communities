import sys

# Check if the file path is provided as an argument
if len(sys.argv) < 3:
    print("Please provide a file path as an argument")
    sys.exit()

# Get the file path from the command line arguments
file_path = sys.argv[1]
d = {}
# Open the file for reading
fout = open(sys.argv[2], 'w')
with open(file_path, "r") as f:
    idx = 0
    # Read the file line by line
    for line in f:
        if line.startswith("#"):
            continue
        endpoints = line.strip().split('\t')
        if endpoints[0] not in d:
            d[endpoints[0]] = idx
            idx += 1
        if endpoints[1] not in d:
            d[endpoints[1]] = idx
            idx += 1
        # Replace all tabs with whitespace
        line = line.replace("\t", " ")
        # Print the modified line
        fout.write(str(d[endpoints[0]])+" "+str(d[endpoints[1]])+'\n')
        
file_path = sys.argv[3]
fout = open(sys.argv[4], 'w')
size_thre = int(sys.argv[5]) if len(sys.argv)>5 else 2
with open(file_path, "r") as f:
    # Read the file line by line
    for line in f:
        if line.startswith("#"):
            continue
        nodes = line.strip().split('\t')
        if len(nodes)<size_thre:
            continue
        # Print the modified line
        fout.write("\t".join([str(d[nodes[i]]) for i in range(len(nodes))])+'\n')
print(idx)
fout.close()