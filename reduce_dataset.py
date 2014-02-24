import sys

def main():
    ids = set([line.strip() for line in open("data/annotation.dat")])

    out_file = open("data/disc_set/discovery_ExpressionMatrix_red.txt", "wb")

    with open("data/disc_set/discovery_ExpressionMatrix.txt", "r") as f:
        header = f.readline()
        out_file.write(header)
        for record in f:
            rec_id = record.split()[0]
            if rec_id in ids:
                out_file.write(record)

    out_file.close()

main()
