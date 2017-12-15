#!/usr/bin/python
import argparse
import sys
import math

class PhysioCal():
    def __init__(self):
        self.data = []
        
    def fasta_reader(self,inputpath):
         fileinput = open(inputpath, "r")
         line = fileinput.readline()
         header = ''
         sequence = ''
         while line:
              line = line.rstrip('\n')
              if '>' in line:
                   header = line
              else:
                   sequence = sequence + line
              line = fileinput.readline()
         print(line)
         print()
         return sequence
         fileinput.close()

    def nucleotide(self,sequence):
         seq=str.upper(sequence)
         length = len(seq)
         A = seq.count('A')
         T = seq.count('T')
         G = seq.count('G')
         C = seq.count('C')
         extinction_C=((A*15.4)+(C*7.2)+(G*11.5)+(T*8.7))*0.9*1000
         m_weight = ((A * 313.21) + (T * 304.2) + (G * 329.21) + (C * 289.18)) - 61.96
         if length < 14:
              temp1 = 4 * (G + C) + 2 * (A + T)
              print('Theoratical Basic temp. is:%.1f C' % temp1)
         elif length >= 14:
              temp2 = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)
              print('Theoratical temp. is:%.1f C' % temp2)

         print('Lenngth of sequence is %d bP:' % length)
         print('Anhydrous Molecular weight:%.1f daltons' % m_weight)
         print('No. of Gaunine in sequence are:%d' % G)
         print('No. of Cytocine in sequence are:%d' % C)
         print('No. of Adnine in sequence are:%d' % A)
         print('No. of Thymine in sequence are:%d' % T)
         print('GC content of seq. is:%d' % ((G + C) / length * 100), '%')
         print('Extenction Coefficient of sequence is: %d /M /cm' %extinction_C)

         '''variable for each IUPAC aminoacid codon'''
         f1 = 0;f2 = 0;l1 = 0;l2 = 0;l3 = 0;l4 = 0;l5 = 0;l6 = 0;s1 = 0;s2 = 0;s3 = 0
         s4 = 0;s5 = 0;s6 = 0;r1 = 0;r2 = 0;r3 = 0;r4 = 0;r5 = 0;r6 = 0;p1 = 0;p2 = 0
         p3 = 0;p4 = 0;v1 = 0;v2 = 0;v3 = 0;v4 = 0;t1 = 0;t2 = 0;t3 = 0;t4 = 0;a1 = 0
         a2 = 0;a3 = 0;a4 = 0;g1 = 0;g2 = 0;g3 = 0;g4 = 0;i1 = 0;i2 = 0;i3 = 0
         stop1 = 0;stop2 = 0;stop3 = 0;y1 = 0;y2 = 0;h1 = 0;h2 = 0;q1 = 0;q2 = 0;n1 = 0
         n2 = 0;k1 = 0;k2 = 0;d1 = 0;d2 = 0;e1 = 0;e2 = 0;c1 = 0;c2 = 0;w1 = 0;m = 0

         x=0
         size=len(seq)
         while (x < size):
              if ((x + 3) <= size):
                   if (seq[x:x + 3] == "TTT"):
                        f1 = f1 + 1
                   elif (seq[x:x + 3] == "TTC"):
                        f2 = f2 + 1
                   elif (seq[x:x + 3] == "TTA"):
                        l1 = l1 + 1
                   elif (seq[x:x + 3] == "TTG"):
                        l2 = l2 + 1
                   elif (seq[x:x + 3] == "CTT"):
                        l3 = l3 + 1
                   elif (seq[x:x + 3] == "CTC"):
                        l4 = l4 + 1
                   elif (seq[x:x + 3] == "CTA"):
                        l5 = l5 + 1
                   elif (seq[x:x + 3] == "CTG"):
                        l6 = l6 + 1
                   elif (seq[x:x + 3] == "ATT"):
                        i1 = i1 + 1
                   elif (seq[x:x + 3] == "ATC"):
                        i2 = i2 + 1
                   elif (seq[x:x + 3] == "ATA"):
                        i3 = i3 + 1
                   elif (seq[x:x + 3] == "ATG"):
                        m = m + 1
                   elif (seq[x:x + 3] == "GTT"):
                        v1 = v1 + 1
                   elif (seq[x:x + 3] == "GTC"):
                        v2 = v2 + 1
                   elif (seq[x:x + 3] == "GTA"):
                        v3 = v3 + 1
                   elif (seq[x:x + 3] == "GTG"):
                        v4 = v4 + 1
                   elif (seq[x:x + 3] == "TCT"):
                        s1 = s1 + 1
                   elif (seq[x:x + 3] == "TCC"):
                        s2 = s2 + 1
                   elif (seq[x:x + 3] == "TCA"):
                        s3 = s3 + 1
                   elif (seq[x:x + 3] == "TCG"):
                        s4 = s4 + 1
                   elif (seq[x:x + 3] == "CCT"):
                        p1 = p1 + 1
                   elif (seq[x:x + 3] == "CCC"):
                        p2 = p2 + 1
                   elif (seq[x:x + 3] == "CCA"):
                        p3 = p3 + 1
                   elif (seq[x:x + 3] == "CCG"):
                        p4 = p4 + 1
                   elif (seq[x:x + 3] == "ACT"):
                        t1 = t1 + 1
                   elif (seq[x:x + 3] == "ACC"):
                        t2 = t2 + 1
                   elif (seq[x:x + 3] == "ACA"):
                        t3 = t3 + 1
                   elif (seq[x:x + 3] == "ACG"):
                        t4 = t4 + 1
                   elif (seq[x:x + 3] == "GCT"):
                        a1 = a1 + 1
                   elif (seq[x:x + 3] == "GCC"):
                        a2 = a2 + 1
                   elif (seq[x:x + 3] == "GCA"):
                        a3 = a3 + 1
                   elif (seq[x:x + 3] == "GCG"):
                        a4 += 1
                   elif (seq[x:x + 3] == "TAT"):
                        y1 = y1 + 1
                   elif (seq[x:x + 3] == "TAC"):
                        y2 = y2 + 1
                   elif (seq[x:x + 3] == "TAA"):
                        stop1 = stop1 + 1
                   elif (seq[x:x + 3] == "TAG"):
                        stop2 = stop2 + 1
                   elif (seq[x:x + 3] == "CAT"):
                        h1 = h1 + 1
                   elif (seq[x:x + 3] == "CAC"):
                        h2 = h2 + 1
                   elif (seq[x:x + 3] == "CAA"):
                        q1 = q1 + 1
                   elif (seq[x:x + 3] == "CAG"):
                        q2 = q2 + 1
                   elif (seq[x:x + 3] == "AAT"):
                        n1 = n1 + 1
                   elif (seq[x:x + 3] == "AAC"):
                        n2 = n2 + 1
                   elif (seq[x:x + 3] == "AAA"):
                        k1 = k1 + 1
                   elif (seq[x:x + 3] == "AAG"):
                        k2 = k2 + 1
                   elif (seq[x:x + 3] == "GAT"):
                        d1 = d1 + 1
                   elif (seq[x:x + 3] == "GAC"):
                        d2 = d2 + 1
                   elif (seq[x:x + 3] == "GAA"):
                        e1 = e1 + 1
                   elif (seq[x:x + 3] == "GAG"):
                        e2 = e2 + 1
                   elif (seq[x:x + 3] == "TGT"):
                        c1 = c1 + 1
                   elif (seq[x:x + 3] == "TGC"):
                        c2 = c2 + 1
                   elif (seq[x:x + 3] == "TGA"):
                        stop3 = stop3 + 1
                   elif (seq[x:x + 3] == "TGG"):
                        w1 = w1 + 1
                   elif (seq[x:x + 3] == "CGT"):
                        r1 = r1 + 1
                   elif (seq[x:x + 3] == "CGC"):
                        r2 = r2 + 1
                   elif (seq[x:x + 3] == "CGA"):
                        r3 = r3 + 1
                   elif (seq[x:x + 3] == "CGG"):
                        r4 = r4 + 1
                   elif (seq[x:x + 3] == "AGT"):
                        s5 = s5 + 1
                   elif (seq[x:x + 3] == "AGC"):
                        s6 = s6 + 1
                   elif (seq[x:x + 3] == "AGA"):
                        r5 = r5 + 1
                   elif (seq[x:x + 3] == "AGG"):
                        r6 = r6 + 1
                   elif (seq[x:x + 3] == "GGT"):
                        g1 = g1 + 1
                   elif (seq[x:x + 3] == "GGC"):
                        g2 = g2 + 1
                   elif (seq[x:x + 3] == "GGA"):
                        g3 = g3 + 1
                   elif (seq[x:x + 3] == "GGG"):
                        g4 = g4 + 1
              x = x + 3
         print('Genome Biasness')
         print("Ala \t GCG \t %d" % a4)
         print("Ala \t GCA \t %d" % a3)
         print("Ala \t GCT \t %d" % a1)
         print("Ala \t GCC \t %d" % a2)
         print("Cys \t TGT \t %d" % c1)
         print("Cys \t TGC \t %d" % c2)
         print("Asp \t GAT \t %d" % d1)
         print("Asp \t GAC \t %d" % d2)
         print("Glu \t GAG \t %d" % e2)
         print("Glu \t GAA \t %d" % e1)
         print("Phe \t TTT \t %d" % f1)
         print("Phe \t TTC \t %d" % f2)
         print("Gly \t GGG \t %d" % g4)
         print("Gly \t GGA \t %d" % g3)
         print("Gly \t GGT \t %d" % g1)
         print("Gly \t GGC \t %d" % g2)
         print("His \t CAT \t %d" % h1)
         print("His \t CAC \t %d" % h2)
         print("Ile \t ATA \t %d" % i3)
         print("Ile \t ATT \t %d" % i1)
         print("Ile \t ATC \t %d" % i2)
         print("Lys \t AAG \t %d" % k2)
         print("Lys \t AAA \t %d" % k1)
         print("Leu \t TTG \t %d" % l2)
         print("Leu \t TTA \t %d" % l1)
         print("Leu \t CTG \t %d" % l6)
         print("Leu \t CTA \t %d" % l5)
         print("Leu \t CTT \t %d" % l3)
         print("Leu \t CTC \t %d" % l4)
         print("Met \t ATG \t %d" % m)
         print("Asn \t AAT \t %d" % n1)
         print("Asn \t AAC \t %d" % n2)
         print("Pro \t CCG \t %d" % p4)
         print("Pro \t CCA \t %d" % p3)
         print("Pro \t CCT \t %d" % p1)
         print("Pro \t CCC \t %d" % p2)
         print("Glu \t CAG \t %d" % q2)
         print("Glu \t CAA \t %d" % q1)
         print("Arg \t AGG \t %d" % r6)
         print("Arg \t AGA \t %d" % r5)
         print("Arg \t CGG \t %d" % r4)
         print("Arg \t CGA \t %d" % r3)
         print("Arg \t CGT \t %d" % r1)
         print("Arg \t CGC \t %d" % r2)
         print("Ser \t AGT \t %d" % s5)
         print("Ser \t AGC \t %d" % s6)
         print("Ser \t TCG \t %d" % s4)
         print("Ser \t TCA \t %d" % s3)
         print("Ser \t TCT \t %d" % s1)
         print("Ser \t TCC \t %d" % s2)
         print("Thr \t ACG \t %d" % t4)
         print("Thr \t ACA \t %d" % t3)
         print("Thr \t ACT \t %d" % t1)
         print("Thr \t ACC \t %d" % t2)
         print("Val \t GTG \t %d" % v4)
         print("Val \t GTA \t %d" % v3)
         print("Val \t GTT \t %d" % v1)
         print("Val \t GTC \t %d" % v2)
         print("Try \t TGG \t %d" % w1)
         print("Tyr \t TAT \t %d" % y1)
         print("Tyr \t TAC \t %d" % y2)
         print("End \t TGA \t %d" % stop3)
         print("End \t TAG \t %d" % stop2)
         print("End \t TAA \t %d" % stop1)

    #protein calculation function
    def amino_acid(self,sequence):
         length = len(sequence)
         '''count individual amino acid'''
         G = sequence.count('G')  # Glysine
         A = sequence.count('A')  # Alanine
         L = sequence.count('L')  # Leucine
         M = sequence.count('M')  # Methionine
         F = sequence.count('F')  # Phenylalanine
         W = sequence.count('W')  # Trypthphan
         K = sequence.count('K')  # Lysine
         Q = sequence.count('Q')  # Glutamine
         E = sequence.count('E')  # Glutamic acid
         S = sequence.count('S')  # Serine
         P = sequence.count('P')  # Proline
         V = sequence.count('V')  # Valine
         I = sequence.count('I')  # Isoleusine
         C = sequence.count('C')  # Cystine
         Y = sequence.count('Y')  # Tyrosine
         H = sequence.count('H')  # Histidine
         R = sequence.count('R')  # Argnine
         N = sequence.count('N')  # Asparagine
         D = sequence.count('D')  # Aspartic Acid
         T = sequence.count('T')  # Threonine

         '''Hydropthy values calculation'''
         a = A * 1.800  # Alanine
         l = L * 3.800  # Leucine
         m = M * 1.900  # Methionine
         f = F * 2.800  # Phenylalanine
         v = V * 4.200  # Valine
         i = I * 4.500  # Isoleusine
         c = C * 2.500  # Cystine
         g = G * (-0.400)  # Glysine
         w = W * (-0.900)  # Trypthphan
         k = K * (-3.900)  # Lysine
         q = Q * (-3.500)  # Glutamine
         e = E * (-3.500)  # Glutamic acid
         s = S * (-0.800)  # Serine
         p = P * (-1.600)  # Proline
         y = Y * (-1.300)  # Tyrosine
         h = H * (-3.200)  # Histidine
         r = R * (-4.500)  # Argnine
         n = N * (-3.500)  # Asparagine
         d = D * (-3.500)  # Aspartic Acid
         t = T * (-0.700)  # Threonine
         total = a + l + m + f + v + i + c + g + w + k + q + e + s + p + y + h + r + n + d + t

         '''molecular weight calculation'''
         w_a = A * 89.0935   # Alanine
         w_l = L * 131.1736  # Leucine
         w_m = M * 149.2124  # Methionine
         w_f = F * 165.1900  # Phenylalanine
         w_v = V * 117.1469  # Valine
         w_i = I * 131.1736  # Isoleusine
         w_c = C * 121.1590  # Cysteine
         w_g = G * 75.0669   # Glysine
         w_w = W * 204.2262  # Trypthphan
         w_k = K * 146.1882  # Lysine
         w_q = Q * 146.1451  # Glutamine
         w_e = E * 147.1299  # Glutamic acid
         w_s = S * 105.0930  # Serine
         w_p = P * 115.1310  # Proline
         w_y = Y * 181.1894  # Tyrosine
         w_h = H * 155.1552  # Histidine
         w_r = R * 174.2017  # Argnine
         w_n = N * 132.1184  # Asparagine
         w_d = D * 133.1032  # Aspartic Acid
         w_t = T * 119.1197  # Threonine
         total2 = w_a + w_l + w_m + w_f + w_v + w_i + w_c + w_g + w_w + w_k + w_q + w_e + w_s + w_p + w_y + w_h + w_r + w_n + w_d + w_t
         h2o = 18.015
         molecular_weight = total2 - h2o * (length - 1)

         '''Protein half life calculation'''
         # mammalian #Yeast #E.coli
         hlife = {
           'A': ['4.4 hour', '>20hour', '>10 hour'],
           'R': ['1 hour', '2 min', '2  min'],
           'N': ['1.4 hour', '3 min', '>10 hour'],
           'D': ['1.1 hour', '3 min', '>10 hour'],
           'C': ['1.2 hour', '>20 hour', '>10 hour'],
           'Q': ['0.8 hour', '10 min', '>10 hour'],
           'E': ['1 hour', '30 min', '>10 hour'],
           'G': ['30 hour', '>20 hour', '>10 hour'],
           'H': ['3.5 hour', '10 min', '>10 hour'],
           'I': ['20 hour', '30 min', '>10 hour'],
           'L': ['5.5 hour', '3 min', '2 min'],
           'K': ['1.3 hour', '3 min', '2 min'],
           'M': ['30 hour', '>20 hour', '>10 hour'],
           'F': ['1.1 hour', '3 min', '2 min'],
           'P': ['>20 hour', '>20 hour', '-'],
           'S': ['1.9 hour', '>20 hour', '>10 hour'],
           'T': ['7.2 hour', '>20 hour', '>10 hour'],
           'W': ['2.8 hour', '3 min', '2 min'],
           'Y': ['2.8 hour', '10 min', '2 min'],
           'V': ['100 hour', '>20 hour', '>10 hour'],
           }
         '''Isoelectric point'''
         Asp = 'D'
         Glu = 'E'
         Cys = 'C'
         Tyr = 'Y'
         His = 'H'
         Lys = 'K'
         Arg = 'R'
         NQ, pH = 0.0, 0.0
         QN1, QN2, QN3, QN4, QN5, QP1, QP2, QP3, QP4 = 0, 0, 0, 0, 0, 0, 0, 0, 0

         var = 1
         while 1:
           QN1 = -1 / (1.0 + pow(10, (3.65 - pH)))
           QN2 = -sequence.count("D") / (1.0 + pow(10, (3.9 - pH)))
           QN3 = -sequence.count("E") / (1.0 + pow(10, (4.07 - pH)))
           QN4 = -sequence.count("C") / (1.0 + pow(10, (8.18 - pH)))
           QN5 =-sequence.count("Y") / (1.0 + pow(10, (10.46 - pH)))
           QP1 = sequence.count("H") / (1.0 + pow(10, (pH - 6.04)))
           QP2 = 1 / (1.0 + pow(10, (pH - 8.2)))
           QP3 = sequence.count("K") / (1.0 + pow(10, (pH - 10.54)))
           QP4 = sequence.count("R") / (1.0 + pow(10, (pH - 12.48)))

           NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4

           if (pH >= 14.0):
               # you should never see this message
               print("Something is wrong - pH is higher than 14")
               break
           elif (NQ <= 0):
               # if this condition is true we can stop calculate
               break
           pH += 0.01  # if not increase pH
         '''extinction coefficient'''
         ext_Trp = 5500
         ext_Tyr = 1490
         ext_Cys = 125
         numb_Tyr = Y
         numb_Trp = W
         numb_Cys = C / 2
         ext_coef = (numb_Tyr * ext_Tyr) + (numb_Trp * ext_Trp) + (int(numb_Cys) * ext_Cys)
         absorb = ext_coef / molecular_weight
         ext_coef2 = (numb_Tyr * ext_Tyr) + (numb_Trp * ext_Trp)
         absorb2 = ext_coef2 / molecular_weight

         '''results display'''
         print("\t\tAA_composition")
         print('Aminoacid \t total number \t %age')
         print("Gly (G) \t  %d \t\t %.2f" % (G, (G / length) * 100))
         print("Ala (A) \t  %d \t\t %.2f" % (A, (A / length) * 100))
         print("Leu (L) \t  %d \t\t %.2f" % (L, (L / length) * 100))
         print("Met (M) \t  %d \t\t %.2f" % (M, (M / length) * 100))
         print("Phe (F) \t  %d \t\t %.2f" % (F, (F / length) * 100))
         print("Trp (W) \t  %d \t\t %.2f" % (W, (W / length) * 100))
         print("Lys (K) \t  %d \t\t %.2f" % (K, (K / length) * 100))
         print("Gln (Q) \t  %d \t\t %.2f" % (Q, (Q / length) * 100))
         print("Glu (E) \t  %d \t\t %.2f" % (E, (E / length) * 100))
         print("Ser (S) \t  %d \t\t %.2f" % (S, (S / length) * 100))
         print("Pro (P) \t  %d \t\t %.2f" % (P, (P / length) * 100))
         print("Val (V) \t  %d \t\t %.2f" % (V, (V / length) * 100))
         print("Ile (I) \t  %d \t\t %.2f" % (I, (I / length) * 100))
         print("Cys (C) \t  %d \t\t %.2f" % (C, (C / length) * 100))
         print("Tyr (Y) \t  %d \t\t %.2f" % (Y, (Y / length) * 100))
         print("His (H) \t  %d \t\t %.2f" % (H, (H / length) * 100))
         print("Arg (R) \t  %d \t\t %.2f" % (R, (R / length) * 100))
         print("Asn (N) \t  %d \t\t %.2f" % (N, (N / length) * 100))
         print("Asp (D) \t  %d \t\t %.2f" % (D, (D / length) * 100))
         print("Thr (T) \t  %d \t\t %.2f" % (T, (T / length) * 100))
         print()
         print("Protein Length is:%d" % length)
         print("Total negetively charged residues (Asp + Glu): %d" % (D + E))
         print("Total positively charged residues (Arg + Lys): %d" % (R + K))
         print("Histidine is considered Nutral because pI=7.37")
         print("Total polar uncharged residues: %d" % (S + T + C + P + N + Q))
         print("Total nonPolar Aliphatic Residues: %d" % (G + A + V + L + I + M))
         print("Total Aromatic residues: %d" % (F + Y + W))
         print("Grand Average of Hydrophobicity (GRAVY):%.3f " % (total / length))
         print("Molecular Weight of protein: %.2f kDa" % (molecular_weight / 1000))

         if (sequence[0] == 'A'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['A'])
         elif (sequence[0] == 'R'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['R'])
         elif (sequence[0] == 'N'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['N'])
         elif (sequence[0] == 'D'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['D'])
         elif (sequence[0] == 'C'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['C'])
         elif (sequence[0] == 'Q'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['Q'])
         elif (sequence[0] == 'E'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['E'])
         elif (sequence[0] == 'G'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['G'])
         elif (sequence[0] == 'H'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['H'])
         elif (sequence[0] == 'I'):
               print('protein half life in Mammalian Yeast and E.coli is:\n', hlife['I'])
         elif (sequence[0] == 'L'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['L'])
         elif (sequence[0] == 'K'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['K'])
         elif (sequence[0] == 'M'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['M'])
         elif (sequence[0] == 'F'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['F'])
         elif (sequence[0] == 'P'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['P'])
         elif (sequence[0] == 'S'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['S'])
         elif (sequence[0] == 'T'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['T'])
         elif (sequence[0] == 'W'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['W'])
         elif (sequence[0] == 'Y'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['Y'])
         elif (sequence[0] == 'V'):
               print('protein half life in Mammalian Yeast and E.coli is:\n ', hlife['V'])
         print("Protein isoelectric point=%.2f" % pH)
         print('Extinction coefficients are in units of  M-1 cm-1, at 280 nm measured in water')
         print('assuming all pairs of Cys residues form cystines')
         print('\tExt. coefficient: %d' % ext_coef)
         print('\tAbs 0.1% (=1g/l):', '%.3f' % absorb)
         print('assuming all pairs of Cys residues are reduced')
         print('\tExt. coefficient: %d' % ext_coef2)
         print('\tAbs 0.1% (=1g/l):', '%.3f' % absorb2)
    def about(self):
        print( '===============================================================================\n' )
        print( '\t\tPhysioCal:A Physical Properties Calculation Script\n' )
        print( 'AUTHOR: \t\tUmer Farooq, UmerBaig7@gmail.com' )
        print( 'COPYRIGHTS: \t\tUmer Farooq\n' )
        print( '\t\t\tDecember 2017' )
        print( '===============================================================================\n' )  

def run(args):
    obj = PhysioCal()
    fin=args.input
    seq=obj.fasta_reader(fin)
    op=args.operation
    obj.about()
    if (op == 'N' or op == 'n'):
        with open (args.output,"w") as fout:
            sys.stdout=fout
            obj.nucleotide(seq)
            fout.close()
    elif (op == 'A' or op == 'a'):
        with open (args.output,"w") as fout:
            sys.stdout=fout
            obj.amino_acid(seq)
            fout.close()

def main():
    parser=argparse.ArgumentParser(description="PhysioCal: Physical Properties Calculation Script")
    parser.add_argument("-in",help="sequence input file", dest="input", type=str, required=True)
    parser.add_argument("-out",help="results output file", dest="output", type=str, required=True)
    parser.add_argument("-op",help="operation to perform,N for nucleotide, A for protein ", dest="operation", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
