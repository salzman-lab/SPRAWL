import pandas as pd
import glob
import os

rename_dict = {
    0: 'PTPRC', 1: 'PECAM1', 2: 'KDR', 3: 'CDH5', 4: 'CD19', 5: 'MS4A1', 6: 'CD79A', 7: 'CD79B', 8: 'BANK1', 9: 'BLK', 10: 'TNFRSF13C',
    11: 'TNFRSF13B', 12: 'PAX5', 13: 'FAIM3', 14: 'S100A2', 15: 'KRT5', 16: 'KRT19', 17: 'SFN', 18: 'MIR205HG', 19: 'S100A14', 20: 'MUC4',
    21: 'SELE', 22: 'CNKSR3', 23: 'ADAMTS9', 24: 'ACKR1', 25: 'GLYAT', 26: 'AK4', 27: 'BBOX1', 28: 'AQP1', 29: 'PDZK1', 30: 'TMEM27',
    31: 'SLC9A3R2', 32: 'CLDN5', 33: 'BTNL9', 34: 'SOX18', 35: 'VWF', 36: 'FLT1', 37: 'CLEC14A', 38: 'GPIHBP1', 39: 'EGFL7', 40: 'TIE1',
    41: 'LDB2', 42: 'PLVAP', 43: 'GAS1', 44: 'CCDC80', 45: 'ISM1', 46: 'IGF1', 47: 'IGFBP6', 48: 'MFAP5', 49: 'SCARA5', 50: 'FBN1',
    51: 'SFRP2', 52: 'ADGRL4', 53: 'F8', 54: 'PTPRB', 55: 'SHANK3', 56: 'EMCN', 57: 'CAVIN2', 58: 'WDR72', 59: 'AQP2', 60: 'KCNE1',
    61: 'SCNN1A', 62: 'SOX17', 63: 'PALMD', 64: 'TINAGL1', 65: 'NPDC1', 66: 'CRIP2', 67: 'HYAL2', 68: 'MMRN2', 69: 'GJA4', 70: 'KCNJ1',
    71: 'NUDT4', 72: 'MAL', 73: 'NOTCH3', 74: 'ACTA2', 75: 'MYL9', 76: 'MYH11', 77: 'LMOD1', 78: 'DES', 79: 'VIM', 80: 'FMO1',
    81: 'ALOX5AP', 82: 'C1QC', 83: 'C1QA', 84: 'C1QB', 85: 'MRC1', 86: 'TREM1', 87: 'TREM2', 88: 'KLRA2', 89: 'BCL11B', 90: 'CD37',
    91: 'CD247', 92: 'KLRB1F', 93: 'KLRB1C', 94: 'GZMB', 95: 'DERL3', 96: 'JCHAIN', 97: 'WT1', 98: 'ACTN4', 99: 'SYNPO', 100: 'DAG1',
    101: 'FOXC1', 102: 'PODXL', 103: 'MME', 104: 'CD3E', 105: 'IL7R', 106: 'CD3G', 107: 'LCK', 108: 'CD4', 109: 'CD8A', 110: 'SPIB',
    111: 'BCL11A', 112: 'CD22', 113: 'TRAT1', 114: 'CD6', 115: 'GZMK', 116: 'CD8B1', 117: 'DNASE1L3', 118: 'ENG', 119: 'STAB2', 120: 'GPR182',
    121: 'SPARC', 122: 'EHD3', 123: 'NRP1', 124: 'FCGR2B', 125: 'STAB1', 126: 'IGFBP4', 127: 'OIT3', 128: 'FLT4', 129: 'GPR116', 130: 'MAF',
    131: 'ESAM', 132: 'COL1A2', 133: 'COLEC11', 134: 'DCN', 135: 'SOD3', 136: 'PRELP', 137: 'PTH1R', 138: 'TCF21', 139: 'RELN', 140: 'ECM1',
    141: 'BGN', 142: 'HMGCS2', 143: 'VTN', 144: 'KNG1', 145: 'UOX', 146: 'RDH7', 147: 'CES3A', 148: 'CYP2F2', 149: 'PCK1', 150: 'HAL',
    151: 'CDH1', 152: 'GLUL', 153: 'GULO', 154: 'OAT', 155: 'MARCO', 156: 'CCR2', 157: 'CX3CR1', 158: 'EMR1', 159: 'CLEC4F', 160: 'CD68',
    161: 'IRF7', 162: 'CFP', 163: 'CTSS', 164: 'CD5L', 165: 'FOLR2', 166: 'SDC3', 167: 'VSIG4', 168: 'CSF1R', 169: 'SLC40A1', 170: 'CCL6',
    171: 'CD3D', 172: 'KLRK1', 173: 'NAPSA', 174: 'CORO1A', 175: 'LGALS3', 176: 'PLAC8', 177: 'LSP1', 178: 'MXD1', 179: 'MMP9', 180: 'CSF3R',
    181: 'HDC', 182: 'CXCR2', 183: 'IL1B', 184: 'PGLYRP1', 185: 'CD44', 186: 'Arg2', 187: 'Elane', 188: 'ZAP70', 189: 'CXCR6', 190: 'PRF1',
    191: 'CTSW', 192: 'IL2RB', 193: 'TXK', 194: 'GIMAP4', 195: 'LY6C2', 196: 'ATP1B3', 197: 'KLF2', 198: 'MS4A6B', 199: 'SIGLECH', 200: 'CYBASC3',
    201: 'RNASE6', 202: 'IRF8', 203: 'DNAJC7', 204: 'TSPAN13', 205: 'TCF4', 206: 'TLR7', 207: 'ITGAX', 208: 'EPCAM', 209: 'KRT7', 210: 'RAMP3',
    211: 'TM4SF1', 212: 'FABP4', 213: 'COL4A2', 214: 'CAV1', 215: 'CD300LG', 216: 'COL4A1', 217: 'PLTP', 218: 'PTRF', 219: 'SPARCL1', 220: 'CD34',
    221: 'APOLD1', 222: 'EPAS1', 223: 'GNAI2', 224: 'CD81', 225: 'CALM1', 226: 'H2-AA', 227: 'H2-AB1', 228: 'LAPTM5', 229: 'ARHGDIB', 230: 'PTPN6',
    231: 'LY86', 232: 'CYBB', 233: 'PLD4', 234: 'CD53', 235: 'SELPLG', 236: 'Scg3', 237: 'Scg5', 238: 'GPX3', 239: 'SFRP5', 240: 'PCSK2',
    241: 'DBPHT2', 242: 'PAPPA2', 243: 'CPE', 244: 'PCSK1N', 245: 'APLP1', 246: 'CHGA', 247: 'PEG3', 248: 'SCG2', 249: 'G6PC2', 250: 'ABCC8',
    251: 'SLC30A8', 252: 'PRLR', 253: 'PTPRN2', 254: 'PPP1R1A', 255: 'ATP2A2', 256: 'PRSS53', 257: 'ERO1LB', 258: 'FAM151A', 259: 'SLC2A2', 260: 'NKX6-1',
    261: 'PDX1', 262: 'MAFA', 263: 'PNLIP', 264: 'CPB1', 265: 'CEL', 266: 'CPA1', 267: '2210010C04RIK', 268: 'PNLIPRP1', 269: 'CELA1', 270: 'PNLIPRP2',
    271: 'HHEX', 272: 'NEUROG3', 273: 'CHGB', 274: 'FAM159B', 275: 'TM4SF4', 276: 'CLDN3', 277: 'DCDC2A', 278: '1700011H14RIK', 279: 'DSG2', 280: 'KRT18',
    281: 'CLDN7', 282: 'SPP1', 283: 'KRT8', 284: 'CLU', 285: 'ATP1B1', 286: 'HNF1B', 287: 'KCNK16', 288: 'TSPAN8', 289: 'BAMBI', 290: 'ACE2',
    291: 'SERPING1', 292: 'MFAP4', 293: 'MXRA8', 294: 'HSD11B1', 295: 'PTGIS', 296: 'PDGFRA', 297: 'LOXL1', 298: 'OGN', 299: 'PDGFRB', 300: 'PCOLCE',
    301: 'GSN', 302: 'LUM', 303: 'SERPINF1', 304: 'CRISPLD2', 305: 'HTRA3', 306: 'AEBP1', 307: 'Blank-0', 308: 'Blank-1', 309: 'Blank-2', 310: 'Blank-3',
    311: 'Blank-4', 312: 'Blank-5', 313: 'Blank-6', 314: 'Blank-7', 315: 'Blank-8', 316: 'Blank-9', 317: 'Blank-10', 318: 'Blank-11', 319: 'Blank-12', 320: 'Blank-13',
    321: 'Blank-14', 322: 'Blank-15', 323: 'Blank-16', 324: 'Blank-17', 325: 'Blank-18', 326: 'Blank-19', 327: 'Blank-20', 328: 'Blank-21', 329: 'Blank-22',
    330: 'Blank-23', 331: 'Blank-24', 332: 'Blank-25', 333: 'Blank-26', 334: 'Blank-27', 335: 'Blank-28', 336: 'Blank-29', 337: 'Blank-30', 338: 'Blank-31',
    339: 'Blank-32', 340: 'Blank-33', 341: 'Blank-34', 342: 'Blank-35', 343: 'Blank-36', 344: 'Blank-37', 345: 'Blank-38', 346: 'Blank-39', 347: 'Blank-40',
    348: 'Blank-41', 349: 'Blank-42', 350: 'Blank-43', 351: 'Blank-44', 352: 'Blank-45', 353: 'Blank-46', 354: 'Blank-47', 355: 'Blank-48', 356: 'Blank-49',
    357: 'Blank-50', 358: 'Blank-51', 359: 'Blank-52', 360: 'Blank-53', 361: 'Blank-54', 362: 'Blank-55', 363: 'Blank-56', 364: 'Blank-57', 365: 'Blank-58',
    366: 'Blank-59', 367: 'Blank-60', 368: 'Blank-61', 369: 'Blank-62', 370: 'Blank-63', 371: 'Blank-64', 372: 'Blank-65', 373: 'Blank-66', 374: 'Blank-67',
    375: 'Blank-68', 376: 'Blank-69', 377: 'Blank-70', 378: 'Blank-71', 379: 'Blank-72', 380: 'Blank-73', 381: 'Blank-74', 382: 'Blank-75', 383: 'Blank-76',
    384: 'Blank-77',
}

csvs = glob.glob('gene_cell/CZB*.csv')
csvs.extend(glob.glob('gene_ont/CZB*.csv'))

for p in csvs:
    df = pd.read_csv(p)
    df['gene'] = df['gene'].map(rename_dict)
    df['gene'] = df['gene'].str.capitalize()
    df.to_csv(p+'.renamed',index=False)



