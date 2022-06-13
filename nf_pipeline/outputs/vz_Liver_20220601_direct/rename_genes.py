import pandas as pd
import glob
import os

rename_dict = {
    '0':'Comt', '1':'Ldha', '2':'Pck1', '3':'Akr1a1', '4':'Ugt2b1', '5':'Acsl5', '6':'Ugt2a3', '7':'Igf1',
    '8':'Errfi1', '9':'Serping1', '10':'Adh4', '11':'Hsd17b2', '12':'Tpi1', '13':'Cyp1a2', '14':'Acsl1', '15':'Akr1d1',
    '16':'Alas1', '17':'Aldh7a1', '18':'G6pc', '19':'Hsd17b12', '20':'Pdhb', '21':'Gpd1', '22':'Cyp7b1', '23':'Pgam1',
    '24':'Hc', '25':'Dld', '26':'Cyp2c23', '27':'Proz', '28':'Acss2', '29':'Psap', '30':'Cald1', '31':'Hsd3b3',
    '32':'Galm', '33':'Cxcl12', '34':'Sardh', '35':'Cebpa', '36':'Aldh3a2', '37':'Gck', '38':'Sdc1', '39':'Pdha1',
    '40':'Npc2', '41':'Hsd17b6', '42':'Aqp1', '43':'Adh7', '44':'Smpdl3a', '45':'Egfr', '46':'Pgm1', '47':'Fasn',
    '48':'Ctsc', '49':'Abcb4', '50':'Fyb', '51':'Alas2', '52':'Gpi1', '53':'Fech', '54':'Lsr', '55':'Psmd3',
    '56':'Gm2a', '57':'Pabpc1', '58':'Cbr4', '59':'Tkt', '60':'Tmem56', '61':'Eif3f', '62':'Cxadr', '63':'Srd5a1',
    '64':'Cyp2c55', '65':'Gnai2', '66':'Gimap6', '67':'Hsd3b2', '68':'Grn', '69':'Rpp14', '70':'Csnk1a1', '71':'Egr1',
    '72':'Mpeg1', '73':'Acsl4', '74':'Hmgb1', '75':'Mpp1', '76':'Lcp1', '77':'Plvap', '78':'Aldh1b1', '79':'Oxsm',
    '80':'Dlat', '81':'Csk', '82':'Mcat', '83':'Hsd17b7', '84':'Epas1', '85':'Eif3a', '86':'Nrp1', '87':'Dek',
    '88':'H2afy', '89':'Bpgm', '90':'Hsd3b6', '91':'Dnase1l3', '92':'Serpinh1', '93':'Tinagl1', '94':'Aldoc', '95':'Cyp2c38',
    '96':'Dpt', '97':'Mrc1', '98':'Minpp1', '99':'Fgf1', '100':'Alcam', '101':'Gimap4', '102':'Cav2', '103':'Eng',
    '104':'Adgre1', '105':'Shisa5', '106':'Csf1r', '107':'Esam', '108':'Unc93b1', '109':'Cnp', '110':'Clec14a', '111':'Kdr',
    '112':'Adpgk', '113':'Gca', '114':'Pkm', '115':'Mkrn1', '116':'Sdc3', '117':'Acaca', '118':'Gpr182', '119':'Bmp2',
    '120':'Tfrc', '121':'Timp3', '122':'Calcrl', '123':'Pfkl', '124':'Wnt2', '125':'Cybb', '126':'Icam1', '127':'Cdh5',
    '128':'Sgms2', '129':'Cd48', '130':'Stk17b', '131':'Tubb6', '132':'Vcam1', '133':'Hgf', '134':'Ramp1', '135':'Arsb',
    '136':'Pld4', '137':'Smarca4', '138':'Fstl1', '139':'Pfkm', '140':'Lhfp', '141':'Lmna', '142':'Cd300lg', '143':'Laptm5',
    '144':'Timp2', '145':'Slc25a37', '146':'Fzd7', '147':'Lyve1', '148':'Acacb', '149':'Cyp1a1', '150':'Eno3', '151':'Cd83',
    '152':'Epcam', '153':'Ltbp4', '154':'Pgm2', '155':'Mertk', '156':'Pth1r', '157':'Itga2b', '158':'Kctd12', '159':'Srd5a3',
    '160':'Bmp5', '161':'Pecam1', '162':'G6pc3', '163':'Cyp17a1', '164':'Stab2', '165':'Cygb', '166':'Col1a2', '167':'Nid1', '168':'Cd44',
    '169':'Ctnnal1', '170':'Ephb4', '171':'Elk3', '172':'Foxq1', '173':'Cxcl14', '174':'Fzd4', '175':'Itgb2', '176':'Tcf7',
    '177':'Srd5a2', '178':'Aldh3b1', '179':'Flt4', '180':'Selp', '181':'Rbpj', '182':'Ep300', '183':'Rhoj', '184':'Fzd1', '185':'Tcf7l2',
    '186':'Ssh2', '187':'Col6a1', '188':'Notch2', '189':'Tcf4', '190':'Tek', '191':'Trim47', '192':'Tent5c', '193':'Ncf1', '194':'Lepr',
    '195':'Pck2', '196':'Lmnb1', '197':'Selplg', '198':'Myh10', '199':'Aldoart1', '200':'Podxl', '201':'Kitl', '202':'Tcf3', '203':'Tspan13',
    '204':'Dll4', '205':'Fzd8', '206':'Lad1', '207':'Procr', '208':'Ccr2', '209':'Akr1c18', '210':'Maml1', '211':'Ms4a1', '212':'Hk3', '213':'Bcam',
    '214':'Fzd5', '215':'Dkk3', '216':'Bank1', '217':'Itgal', '218':'Pgam2', '219':'Axin2', '220':'Pfkp', '221':'Meis2', '222':'Jag1', '223':'Gimap3',
    '224':'Rassf4', '225':'Notch1', '226':'Cd93', '227':'Tet2', '228':'Tcf7l1', '229':'Cd34', '230':'Hvcn1', '231':'Mal', '232':'Itgb7',
    '233':'Wnt4', '234':'Kit', '235':'Gapdhs', '236':'Kcnj16', '237':'Tnfrsf13c', '238':'Hk1', '239':'Pdgfra', '240':'Apobec3', '241':'Slc34a2',
    '242':'Vav1', '243':'Lamp3', '244':'Meis1', '245':'Lck', '246':'Efnb2', '247':'Notch4', '248':'Klrb1c', '249':'Angpt2', '250':'Vwf', '251':'E2f2',
    '252':'Ccr1', '253':'Angpt1', '254':'B4galt6', '255':'Cyp21a1', '256':'Pdpn', '257':'Dll1', '258':'Ammecr1', '259':'Csf3r', '260':'Ndn', '261':'Fgf2',
    '262':'Runx1', '263':'Mpl', '264':'Mecom', '265':'Itgam', '266':'Hoxb4', '267':'Tox', '268':'Prickle2', '269':'Acss1', '270':'Cyp2b9', '271':'Aldh3a1',
    '272':'Bmp7', '273':'Gata2', '274':'Il7r', '275':'Satb1', '276':'Sfrp1', '277':'Eno2', '278':'Mrvi1', '279':'Mki67', '280':'Nes', '281':'Tmod1', '282':'Ace',
    '283':'Gfap', '284':'Tgfb2', '285':'Tomt', '286':'Flt3', '287':'Sult2b1', '288':'Hkdc1', '289':'Notch3', '290':'Cdh11', '291':'Il6', '292':'Hk2', '293':'Mmrn1',
    '294':'Vangl2', '295':'Pou2af1', '296':'Hoxb5', '297':'Jag2', '298':'Aldh3b2', '299':'Gypa', '300':'Lrp2', '301':'Lef1', '302':'Olr1', '303':'Lox', '304':'Txlnb',
    '305':'Slc12a1', '306':'Aldh3b3', '307':'Cxcr2', '308':'Nkd2', '309':'Sult1e1', '310':'Acsl6', '311':'Ddx4', '312':'Ldhc', '313':'Kcnj1', '314':'Acsbg1',
    '315':'Fzd3', '316':'F13a1', '317':'Hsd11b2', '318':'Dkk2', '319':'Hsd17b1', '320':'Fzd2', '321':'Cyp2b23', '322':'Eno4', '323':'Celsr2', '324':'Obscn',
    '325':'Slamf1', '326':'Akap14', '327':'Gnaz', '328':'Cd177', '329':'Tet1', '330':'Cspg4', '331':'Aldoart2', '332':'Cyp2b19', '333':'Ryr2', '334':'Ldhal6b',
    '335':'Acsf3', '336':'Chodl', '337':'Ivl', '338':'Cyp11b1', '339':'Sfrp2', '340':'Dkk1', '341':'Cyp11a1', '342':'1700061G19Rik', '343':'Acsbg2', '344':'Olah',
    '345':'Pdha2', '346':'Hsd17b3', '347':'Blank-0', '348':'Blank-1', '349':'Blank-2', '350':'Blank-3', '351':'Blank-4', '352':'Blank-5', '353':'Blank-6',
    '354':'Blank-7', '355':'Blank-8', '356':'Blank-9', '357':'Blank-10', '358':'Blank-11', '359':'Blank-12', '360':'Blank-13', '361':'Blank-14', '362':'Blank-15',
    '363':'Blank-16', '364':'Blank-17', '365':'Blank-18', '366':'Blank-19', '367':'Blank-20', '368':'Blank-21', '369':'Blank-22', '370':'Blank-23', '371':'Blank-24',
    '372':'Blank-25', '373':'Blank-26', '374':'Blank-27', '375':'Blank-28', '376':'Blank-29', '377':'Blank-30', '378':'Blank-31', '379':'Blank-32', '380':'Blank-33',
    '381':'Blank-34', '382':'Blank-35', '383':'Blank-36', '384':'Blank-37',
}

csvs = glob.glob('gene_cell/*.csv')
csvs.extend(glob.glob('gene_ont/*.csv'))

for p in csvs:
    df = pd.read_csv(p)
    df['gene'] = df['gene'].astype(str).map(rename_dict)
    df.to_csv(p+'.renamed',index=False)



