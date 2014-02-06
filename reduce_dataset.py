

def main():
    ids = set((pid.strip() for pid in open("annotation.dat")))
    out_file = open("disc_set/discovery_ExpressionMatrix_red.txt", "wb")
    
    with open("disc_set/discovery_ExpressionMatrix.txt", "r") as f:
        header = f.readline()
        out_file.write(header)
        for record in f:
            rec_id = record.split()[0]
            if rec_id in ids:
                out_file.write(record)


    out_file.close()


main()
