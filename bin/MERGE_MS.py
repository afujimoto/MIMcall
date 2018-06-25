import sys

#1       890208  890256  (A)n    AAAAAAAAAAACCAAAAACCAAACAAACAAAACAAAACAAAACAAAAC        48;10   49;1    

f = open(sys.argv[1])
MS_dic = {}
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	pos = "\t".join(line_l[0:4])
	MS_dic[pos] = [line_l]

f2 = open(sys.argv[2])
for line in f2:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	pos = "\t".join(line_l[0:4])
	if pos in MS_dic:
		MS_dic[pos].append(line_l)

#1       198849625       198849666       (ATTT)n ATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTA   45;7,   45;3,49;4,

for pos in MS_dic:
	if len(MS_dic[pos]) == 2:
		MS_allele1 = []
		for i in range(5, len(MS_dic[pos][0])):
			MS_allele1.append(MS_dic[pos][0][i])

		MS_allele2 = []
		for i in range(5, len(MS_dic[pos][1])):
			MS_allele2.append(MS_dic[pos][1][i])

		print(pos, MS_dic[pos][0][4], ",".join(MS_allele1), ",".join(MS_allele2), sep="\t")
