#!/usr/bin/env python

## be sure to activate environment where pycairo is installed

import argparse
import re
import cairo


def get_args():
    parser = argparse.ArgumentParser(description="Visualize cassette exon splicing motifs. Outputs an SVG file.")
    parser.add_argument("-f", "--file", help="fasta file with gene sequence(s) to be visualized, exons must be capitalized and introns lowercase", required=True, type=str)
    parser.add_argument("-m", "--motifs", help="text file with motifs one per line, IUPAC compatible", required=True, type=str)
    return parser.parse_args()

# store args
args = get_args()
genes_fasta = args.file
motifs_file = args.motifs

# extract file prefix
spl_name = re.split("\.", genes_fasta)
prefix = spl_name[0]    # store prefix for naming output svg


# constants for figure
HEIGHT_PER_GENE = 150
LEFT_MARGIN = 50
RIGHT_MARGIN = 400      # make big enough to fit legend
FEATURE_HEIGHT = 35     # height of exon and motif rectangles
LABEL_FONT = 14
ALPHA = .75

# rgb colors, will be converted to 0-1 values
colors256 = [
    (58, 52, 235),  # blue-purple
    (50, 168, 82),  # green
    (31, 237, 237), # light blue
    (237, 127, 31), # orange
    (255, 102, 229)   # fuscia
]

# IUPAC encodings
IUPAC = {
    "A":"A",
    "C":"C",
    "G":"G",
    "T":"T",
    "U":"T",  # change U to T
    "W":"[A,T]",
    "S":"[C,G]",
    "M":"[A,C]",
    "K":"[G,T]",
    "R":"[A,G]",
    "Y":"[C,T]",
    "B":"[C,G,T]",
    "D":"[A,G,T]",
    "H":"[A,C,T]",
    "V":"[A,C,G]",
    "N":"[A,C,G,T]",
    "Z":"[]",
}

## Classes -----------------------------------------------------------------------

class FastaHeader:
    '''header object with method to draw itself'''
    def __init__(self, context, gene_num, gene_dict_vals):
    ## Data ##
        self.gnumber = gene_num
        self.context = context
        self.gene_dict_vals = gene_dict_vals

    ## Methods ##
    def draw(self):
        gene_ypos = (self.gnumber*HEIGHT_PER_GENE) - (HEIGHT_PER_GENE/2)
        label_pos = gene_ypos - (HEIGHT_PER_GENE/4)
        self.context.set_source_rgb(0, 0, 0) 
        self.context.set_font_size(LABEL_FONT)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(LEFT_MARGIN, (label_pos))
        self.context.show_text(self.gene_dict_vals[0])
        self.context.stroke()

class Exon:
    '''exon object with method to draw itself'''
    def __init__(self, context, gene_num, gene_dict_vals):
    ## Data ##
        self.gnumber = gene_num
        self.context = context
        self.gene_dict_vals = gene_dict_vals

    ## Methods ##
    def draw(self):
        ypos = self.gnumber * HEIGHT_PER_GENE
        exon_ypos = ypos - (HEIGHT_PER_GENE/2)
        exon_start = LEFT_MARGIN + self.gene_dict_vals[2][0]
        exon_end = LEFT_MARGIN + self.gene_dict_vals[2][0] + self.gene_dict_vals[2][1]
        self.context.set_line_width(FEATURE_HEIGHT)
        self.context.move_to(exon_start, exon_ypos)
        self.context.line_to(exon_end, exon_ypos)
        self.context.stroke()

class Motifs:
    '''motif object with method to draw all motifs in for a gene group'''
    def __init__(self, context, gene_num, gene_dict_vals, colors):
    ## Data ##
        self.gnumber = gene_num
        self.context = context
        self.gene_dict_vals = gene_dict_vals
        self.colors = colors

    def draw(self):
        motif_ypos = (self.gnumber*HEIGHT_PER_GENE) - (HEIGHT_PER_GENE/2)
        all_motifs = []     # list for compiling motif info
        for i in range(len(MOTIF_PATTERNS)):
            mtfs = get_motif_positions(MOTIF_PATTERNS[i][1], MOTIF_PATTERNS[i][0], self.gene_dict_vals[1].upper())    # {IUPAC motif:[(start, end, match),...()]}
            all_motifs.append(mtfs)
            for n in mtfs.values():
                for o in n:     # draw motifs, different motifs have different colors
                    self.context.set_source_rgba(self.colors[i][0], self.colors[i][1], self.colors[i][2], ALPHA)
                    self.context.set_line_width(FEATURE_HEIGHT)
                    self.context.move_to((LEFT_MARGIN+o[0]), motif_ypos)
                    self.context.line_to((LEFT_MARGIN+o[1]), motif_ypos)
                    self.context.stroke()

class GeneGroup:
    '''gene group object with method to draw line with length of the gene and calls draw method for header, exon, and motif objects'''
    def __init__(self, context, gene_num, gene_dict_vals, header, exon, motifs):
    ## Data ##
        self.context = context
        self.gnumber = gene_num
        self.gene_dict_vals = gene_dict_vals
        self.header = header
        self.exon = exon
        self.motifs = motifs
    
    def draw(self):
        gene_ypos = (self.gnumber*HEIGHT_PER_GENE) - (HEIGHT_PER_GENE/2)
        # draw line
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_line_width(1)
        self.context.move_to(LEFT_MARGIN, gene_ypos)
        self.context.line_to((LEFT_MARGIN+self.gene_dict_vals[2][3]), gene_ypos)
        self.context.stroke()

        self.header.draw()
        self.exon.draw()
        self.motifs.draw()


## Functions ---------------------------------------------------------------------

## genes --
def parse_introns_exon(seq):
    '''input sequence with lowercase intron, uppercase exon, lowercase intron,
    returns length of intron1, the exon, and intron2'''
    i1 = 0  # intron
    e = 0   # exon
    i2 = 0   # intron
    total = 0
    for n in seq:
        if e == 0 and n.islower():
            i1+=1
        elif n.isupper():
            e+=1
        elif i1 != 0 and e != 0 and n.islower():
            i2+=1
    total = i1 + e + i2
    return (i1,e,i2, total)

def get_genes(fasta):
    '''parses a fasta file where the introns are lowercase and exons are uppercase,
    returns dictionary where keys are gene1,gene2,etc. and values are [header, sequence, (intron1 len, exon len, intron2 len)]'''
    with open(fasta, "r") as fa:
        genedict = {}       # {gene1: [header, sequence, (intron1 len, exon len, intron2 len)]}
        header = ""
        record = []
        n = 1
        for line in fa:
            line = line.strip()
            if record != []:
                if line[0] != ">":                # if record has something in it and line is a sequence line
                    record.append(line)           # append sequence
                else:                             # if record has something and is a header line, then what is stored in record is a full record
                    header = record[0].strip(">")            # header first element in the record list
                    whole_seq = ""          # concat to get whole sequence
                    for s in record[1:]:
                        whole_seq += s

                    in_ex_lens = parse_introns_exon(whole_seq)    # get introns and exon lengths

                    genedict[n] = [header, whole_seq, in_ex_lens]    # add to dict

                    record = []             # clear record
                    record.append(line)     # add header for the new record
                    n+=1
            else:                    # append first header line
                record.append(line)
            
            # need to do these same steps for the last record
            header = record[0].strip(">")
            whole_seq = ""
            for s in record[1:]:
                whole_seq += s
            in_ex_lens = parse_introns_exon(whole_seq)
            genedict[n] = [header, whole_seq, in_ex_lens]
    return genedict


## motifs --
def get_motif_patterns(txt_file):
    '''input text file where there is 1 IUPAC encoded motif per line,
    converts to uppercase, and returns a list of motif regex patterns'''
    motif_ls = []
    with open(txt_file, "r") as txt:
        n = 1
        for line in txt:
            pattern = ""
            line = line.strip()
            for l in line.upper():      # build regex pattern using IUPAC dict
                pattern += IUPAC[l]
            motif_ls.append((line.upper(), pattern))
            n+=1
    return(motif_ls)

def get_motif_positions(pattern, orig_motif, seq):
    '''input regex pattern generated by get_motif_patterns(), the original IUPAC motif, and sequence to scan,
    converts U's in sequence to T's,
    returns dtionary where keys are the original IUPAC motif and values are the start, end, and seq that matched'''
    motif_spans = []
    if "U" in seq.upper():  # check for U's in sequence
        seq = re.sub("U", "T", seq.upper())
    for x in re.finditer(pattern, seq):
        motif_spans.append((x.start(), x.end(), x.group()))
    return {orig_motif:motif_spans}     # {IUPAC motif:[(start, end, match),...()]}


# figure functions --
def generate_fig_dim(gene_dict):
    '''uses gene information to determine size of the figure,
    returns the width and height'''
    longest_len = 0
    for g in gene_dict.values():    # determine longest gene
        if len(g[1]) > longest_len:
            longest_len = len(g[1])
    fig_width = (longest_len + RIGHT_MARGIN) // 100 * 100    # for legend and round to 100
    fig_height = len(gene_dict) * HEIGHT_PER_GENE   # height depends on number of genes
    return fig_width, fig_height

def color_transf(col_ls):
    '''transforms rgb from 0-256 to 0-1,
    returns list of (r,g,b) tuples as floats'''
    new_col_ls = []
    for c in col_ls:
        new_c = []
        for d in c:
            new_c.append(d/256)
        new_col_ls.append((new_c[0],new_c[1],new_c[2]))
    return new_col_ls

def draw_legend(motif_ls, cols, ctx):
    '''input motif patterns list, colors as floats (0-1), and the context surface,
    draw legend'''
    legend_xstart = width - 300     # choose offset smaller than RIGHT_MARGIN
    legend_ystart = 75
    col_box_size = 20

    # write "Legend"
    ctx.set_source_rgb(0, 0, 0) 
    ctx.set_font_size(col_box_size)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.move_to(legend_xstart, legend_ystart)
    ctx.show_text("Legend")
    ctx.stroke()

    for i in range(len(motif_ls)):
        legend_ypos = legend_ystart + (i*30) + 20
        
        # draw color box
        ctx.set_source_rgba(cols[i][0], cols[i][1], cols[i][2], ALPHA)
        ctx.set_line_width(col_box_size)
        ctx.move_to(legend_xstart, legend_ypos)
        ctx.line_to(legend_xstart+col_box_size, legend_ypos)
        ctx.stroke()

        # write orig motifs
        ctx.set_source_rgb(0, 0, 0) 
        ctx.set_font_size(col_box_size)
        ctx.move_to(legend_xstart+col_box_size+10, legend_ypos+5)
        ctx.show_text(motif_ls[i][0])
        ctx.stroke()


# Main ---------------------------------------------------------------------------

# parse input files
GENES = get_genes(genes_fasta)      # {1:[header, sequence, (intron1 len, exon len, intron2 len)], 2:...}
print(GENES)
MOTIF_PATTERNS = get_motif_patterns(motifs_file)    # [(original motif seq, regex pattern),...]
print(MOTIF_PATTERNS)

# create svg
width, height = generate_fig_dim(GENES)
dec_colors = color_transf(colors256)
surface = cairo.SVGSurface((prefix+"_oop.svg"), width, height)
CONTEXT = cairo.Context(surface)

# store objects
HEADER_OBJS = []
for h in range(len(GENES)):                 # 0 vs 1-based indexing
    HEADER_OBJS.append(FastaHeader(CONTEXT, (h+1), GENES[(h+1)]))
EXON_OBJS = []
for e in range(len(GENES)):
    EXON_OBJS.append(Exon(CONTEXT, (e+1), GENES[(e+1)]))
MOTIFS_OBJS = []
for m in range(len(GENES)):
    MOTIFS_OBJS.append(Motifs(CONTEXT, (m+1), GENES[(m+1)], dec_colors))
GENE_OBJS = []
for g in range(len(GENES)):
    gene_group = GeneGroup(CONTEXT, (g+1), GENES[(g+1)], HEADER_OBJS[g], EXON_OBJS[g], MOTIFS_OBJS[g])
    GENE_OBJS.append(gene_group)

# draw gene groups
for gene_group in GENE_OBJS:
    gene_group.draw()

# draw legend
draw_legend(MOTIF_PATTERNS, dec_colors, CONTEXT)

surface.finish()




