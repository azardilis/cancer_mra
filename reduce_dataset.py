import sys

def read_annotation(annot_file):
    ids = dict()
    for record in open(annot_file):
        pid, tf = record.split(" ")
        ids[pid] = tf.strip()

    return ids

def main():
    ids = read_annotation("annotation.dat")

    out_file = open("data/disc_set/discovery_ExpressionMatrix_red.txt", "wb")
    tf_file = open("data/disc_set/tfs.txt", "wb")
    with open("data/disc_set/discovery_ExpressionMatrix.txt", "r") as f:
        header = f.readline()
        out_file.write(header.strip())
        for record in f:
            rec_id = record.split()[0]
            if rec_id in ids:
                out_file.write(record)
                tf_file.write(ids[rec_id]+"\n")


    out_file.close()
    tf_file.close()

main()
