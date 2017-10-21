import sys
import csv

def go_merge_parser(input_file):

    gene_dic={}
    go_list=[]
    with open(input_file) as p_infile:
        infile = csv.reader(p_infile, delimiter=',')
        gene_tmp = ""
        i = 0
        for rows in infile:
            if "Gene stable ID" in rows:
                continue
            else:
                Gene_ID, Gene_name, GO_id, GO_name, Gene_description, Gene_biotype = rows
                if i == 0:
                    gene_tmp = "\t".join([Gene_ID, Gene_name, Gene_description, Gene_biotype])

                if gene_tmp == "\t".join([Gene_ID, Gene_name, Gene_description, Gene_biotype]):
                    i += 1
                    go = GO_id + "~" + GO_name
                    go_list.append(go)

                else:
                    gene_dic[gene_tmp] = ";".join(go_list)
                    go_list = []

                    go = GO_id + "~" + GO_name
                    go_list.append(go)
                    gene_tmp = "\t".join([Gene_ID, Gene_name, Gene_description, Gene_biotype])
                    i += 1

        gene_dic[gene_tmp] = ";".join(go_list)



    return gene_dic


def main():
    output=open(sys.argv[2],"w")
    dic1=go_merge_parser(sys.argv[1])
    output.write("\t".join(["Gene stable ID","Gene name","Gene description","Gene biotype","GO"])+"\n")
    for gene,go in dic1.items():
        gene_item=gene+"\t"+go+"\n"
        output.write(gene_item)
    output.close()


main()
