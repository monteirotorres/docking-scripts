import os
for folder in os.listdir():
    if os.path.isdir(folder) and folder != "ligands" and folder != "pdbs":
        print(folder)
        for output in os.listdir(folder):
            if output. endswith("txt"):
                with open(os.path.join(folder,output), "r") as fin:
                    for line in fin:
                        pass
                    last_line = line
                    score = last_line.split()[1]
                    with open("Scores.txt", "a") as fout: 
                        fout.write(folder+"\t"+score+"\n")
                    with open("Scores2Merge.txt", "a") as fout:
                        fout.write(folder.split("_")[1]+"\t"+score+"\n")
