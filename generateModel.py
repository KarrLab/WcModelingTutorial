'''
Generates the random components of the example model:
- Random RNA and protein species for each gene in the model
- Transcription, translation, and RNA degradation reactions for each RNA and protein species 
- Biomass production reaction (FBA objective) which is consistent with the specified biomass composition and transcription, translation, and RNA degradation reactions

Author: Jonathan Karr, karr@mssm.edu
Last updated: 3/24/2016
'''

#required libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import random
from openpyxl import Workbook
from openpyxl.styles import Color, PatternFill
import math
import numpy as np

#model parameters
nAvogadro = 6.022e23 #Avogadro constant
cellVol = 4.58e-17 #volume (L)
concAa = 0.5e-3 #[ALA], [ARG], [CYS], ... (M)
concNxp = 1e-3 #[NTP], [NDP], [NMP] (M)
concPPi = 1e-3 #[PPi] (M)
concPi = 5e-3 #[Pi](M)
concH2O = 55 #[H2O] (M)
concH = 1.1220e-08 # [H] (M)
nRnaCopy = 5. #Copy numbers of RNA
nProtCopy = 100. #copy numbers of proteins
CellCycleLength = 10 * 60 * 60. #(s)
RnaHalfLife = 5 * 60. #(s)

#list of to generate random sequences for
genes = [
        {'id': 'RnaPolymerase', 'name': 'RNA polymerase'},
        {'id': 'Ribosome', 'name': 'Ribosome'},
        {'id': 'Rnase', 'name': 'RNase'},
        {'id': 'Adk', 'name': 'Adenylate kinase'},
        {'id': 'Apt', 'name': 'Adenine phosphoribosyltransferase '},
        {'id': 'Cmk', 'name': 'Cytidylate kinase'},
        {'id': 'Eno', 'name': 'Enolase'},
        {'id': 'Fba', 'name': 'Fructose-biphosphate aldolase'},
        {'id': 'Gap', 'name': 'Glyceraldehyde-3-phosphate-dehydrogenase'},
        {'id': 'Pgm', 'name': '2,3-bisphosphoglycerate-independent phosphoglycerate mutase '},
        {'id': 'Gk', 'name': 'Guanylate kinase'},
        {'id': 'Hpt', 'name': 'Hypoxanthine-guanine phosphoribosyltransferase'},
        {'id': 'Ldh', 'name': 'L-lactate dehydrogenase'},
        {'id': 'Nox', 'name': 'Probable NADH oxidase'},
        {'id': 'Pfk', 'name': '6-phosphofructokinase /Phosphohexokinase'},
        {'id': 'PgiB', 'name': 'Glucose-6-phosphate isomerase '},
        {'id': 'Pgk', 'name': 'Phosphoglycerate kinase'},
        {'id': 'Ppa', 'name': 'Inorganic pyrophosphatase '},
        {'id': 'Prs', 'name': 'Ribose-phosphate pyrophosphokinase '},
        {'id': 'Pyk', 'name': 'Pyruvate kinase'},
        {'id': 'Rpe', 'name': 'Probable ribulose-phosphate 3-epimerase'},
        {'id': 'LacA', 'name': 'Probable ribose-5-phosphate isomerase B '},
        {'id': 'Tim', 'name': 'Triosephosphate isomerase '},
        {'id': 'TklB', 'name': 'Transketolase '},
        {'id': 'Udk', 'name': 'Uridine kinase'},
        {'id': 'PyrH', 'name': 'Uridylate kinase '},
        {'id': 'Upp', 'name': 'Uracil phosphoribosyltransferase '},
        {'id': 'Pts', 'name': 'PTS EIIABC (glucose)'},
        {'id': 'PiAbcTransporter', 'name': 'Phosphate ABC transporter'},
        {'id': 'PeptAbcTransporter', 'name': 'Peptide ABC transporter'},
        ]

#1 and 3 letter amino acid codes
aminoAcids = [
    {'id': 'A', 'name': 'ALA', 'charge': 0},
    {'id': 'R', 'name': 'ARG', 'charge': 1},
    {'id': 'D', 'name': 'ASP', 'charge': -1},
    {'id': 'N', 'name': 'ASN', 'charge': 0},
    {'id': 'C', 'name': 'CYS', 'charge': 0},
    {'id': 'Q', 'name': 'GLN', 'charge': 0},
    {'id': 'E', 'name': 'GLU', 'charge': -1},
    {'id': 'G', 'name': 'GLY', 'charge': 0},
    {'id': 'H', 'name': 'HIS', 'charge': 0},
    {'id': 'I', 'name': 'ILE', 'charge': 0},
    {'id': 'L', 'name': 'LEU', 'charge': 0},
    {'id': 'K', 'name': 'LYS', 'charge': 1},
    {'id': 'M', 'name': 'MET', 'charge': 0},
    {'id': 'F', 'name': 'PHE', 'charge': 0},
    {'id': 'P', 'name': 'PRO', 'charge': 0},
    {'id': 'S', 'name': 'SER', 'charge': 0},
    {'id': 'T', 'name': 'THR', 'charge': 0},
    {'id': 'W', 'name': 'TRP', 'charge': 0},
    {'id': 'Y', 'name': 'TYR', 'charge': 0},
    {'id': 'V', 'name': 'VAL', 'charge': 0},
    ]

#1. Generates random RNA and protein species of length 3 * (protLen + 1) and protLen for each gene in the "genes" array. 
#2. Generates transcription, translation, and RNA degradation reactions for each RNA and protein species
#3. Calculates a biomass production reaction (FBA objective) which is consistent with the specified biomass composition (metabolite concentrations) and transcription, translation, and RNA degradation reactions.
#4. Outputs these species and reactions to an Excel file
def run(protLen = 100, startCodon = 'ATG', stopCodon = 'TAG'):
    nGene = len(genes)
    
    '''generate random sequences'''
    dnaSeqs = []
    rnaSeqs = []
    protSeqs = []
    for gene in genes:
        dnaSeq = generateRandomSequence(protLen, startCodon, stopCodon)
        
        #create gene, RNA, and protein sequence records
        dnaSeqs.append(SeqRecord(dnaSeq, id = '%s-Gene' % gene['id'], description = '%s (Gene)' % gene['name']))
        rnaSeqs.append(SeqRecord(dnaSeq.transcribe(), id = '%s-Rna' % gene['id'], description = '%s (RNA)' % gene['name']))       
        protSeqs.append(SeqRecord(Seq(str(dnaSeq.translate(table=4))[0:-1]), id = '%s-Protein' % gene['id'], description = '%s (Protein)' % gene['name']))
        
    '''save sequences and reactions to Excel workbook'''
    wb = Workbook()
    
    specWs = wb.active
    specWs.title = "Species"
    rxnWs = wb.create_sheet(title='Reactions')
        
    #first rows
    specWs['A1'] = 'ID'
    specWs['B1'] = 'Name'
    specWs['C1'] = 'Structure'
    specWs['D1'] = 'Charge'
    specWs['E1'] = 'Type'
    specWs['F1'] = 'Subtype'
    specWs['G1'] = 'Average concentration, cytosol'
    specWs['H1'] = 'Average concentration, extracellular space'
    specWs['I1'] = 'Average count, cytosol'
    specWs['J1'] = 'Average count, extracellular space'
    specWs['K1'] = 'Cross reference source'
    specWs['L1'] = 'Cross reference ID'
    
    rxnWs['A1'] = 'ID'
    rxnWs['B1'] = 'Name'
    rxnWs['C1'] = 'Submodel'
    rxnWs['D1'] = 'Stoichiometry'
    rxnWs['E1'] = 'Enzyme'
    rxnWs['F1'] = 'Rate law'
    rxnWs['G1'] = 'Vmax (1/s, 1/(M*s))'
    rxnWs['H1'] = 'Km (M)'
    rxnWs['I1'] = 'Cross reference source'
    rxnWs['J1'] = 'Cross reference ID'
    
    for iGene, gene in enumerate(genes):
        #RNA species
        specWs['A%d' % (iGene + 2)] = '%s-Rna' % gene['id'] #ID
        specWs['B%d' % (iGene + 2)] = '%s (RNA)' % gene['name'] #name
        specWs['C%d' % (iGene + 2)] = str(rnaSeqs[iGene].seq) #Structure (sequence)
        specWs['D%d' % (iGene + 2)] = -4 * (protLen + 1) * 3 #charge
        specWs['E%d' % (iGene + 2)] = 'RNA' #type
        specWs['F%d' % (iGene + 2)] = '' #subtype
        specWs['G%d' % (iGene + 2)] =  nRnaCopy / cellVol / nAvogadro #cytosol concentration
        specWs['H%d' % (iGene + 2)] =  0 #extracellular concentration
        specWs['I%d' % (iGene + 2)] =  nRnaCopy #cytosol count
        specWs['J%d' % (iGene + 2)] =  0 #extracellular count
        specWs['K%d' % (iGene + 2)] = '' #Cross reference source
        specWs['L%d' % (iGene + 2)] = '' #Cross reference ID
        
        #protein species
        protCharge = 0
        for aa in aminoAcids:
            protCharge += protSeqs[iGene].seq.count(aa['id']) * aa['charge'] - 4 * (2 * len(protSeqs[iGene].seq) + 3 - (len(protSeqs[iGene].seq) - 1))
        
        specWs['A%d' % (iGene + 2 + nGene)] = '%s-Protein' % gene['id'] #ID
        specWs['B%d' % (iGene + 2 + nGene)] = '%s (Protein)' % gene['name'] #name
        specWs['C%d' % (iGene + 2 + nGene)] = str(protSeqs[iGene].seq) #Structure
        specWs['D%d' % (iGene + 2 + nGene)] = protCharge #charge
        specWs['E%d' % (iGene + 2 + nGene)] = 'Protein' #type
        specWs['F%d' % (iGene + 2 + nGene)] = '' #subtype
        specWs['G%d' % (iGene + 2 + nGene)] =  nProtCopy / cellVol / nAvogadro #cytosol concentration
        specWs['H%d' % (iGene + 2 + nGene)] =  0 #extracellular concentration
        specWs['I%d' % (iGene + 2 + nGene)] =  nProtCopy #cytosol count
        specWs['J%d' % (iGene + 2 + nGene)] =  0 #extracellular count
        specWs['K%d' % (iGene + 2 + nGene)] = '' #Cross reference source
        specWs['L%d' % (iGene + 2 + nGene)] = '' #Cross reference ID
        
        #transcription reaction
        nA = rnaSeqs[iGene].seq.count('A')
        
        rxnWs['A%d' % (iGene + 2)] = 'Transcription-%s' % gene['id'] #ID
        rxnWs['B%d' % (iGene + 2)] = 'Transcription (%s)' % gene['name'] #name
        rxnWs['C%d' % (iGene + 2)] = 'Transcription' #submodel
        rxnWs['D%d' % (iGene + 2)] = '[c]: (%d) ATP + (%d) CTP + (%d) GTP + (%d) UTP + H2O ==> %s-Rna + (%d) PPI + H' % ( #stoichometry
            rnaSeqs[iGene].seq.count('A'),
            rnaSeqs[iGene].seq.count('C'),
            rnaSeqs[iGene].seq.count('G'),
            rnaSeqs[iGene].seq.count('U'),
            gene['id'],
            len(rnaSeqs[iGene].seq), )
        rxnWs['E%d' % (iGene + 2)] = 'RnaPolymerase-Protein' #enzyme
        rxnWs['F%d' % (iGene + 2)] = 'Vmax * min(ATP[c], CTP[c], GTP[c], UTP[c]) / (Km + min(ATP[c], CTP[c], GTP[c], UTP[c])) * RnaPol-Protein[c]' #rate law
        rxnWs['G%d' % (iGene + 2)] = math.log(2) * (1 / RnaHalfLife + 1 / CellCycleLength) * 2 * nRnaCopy / nProtCopy #Vmax
        rxnWs['H%d' % (iGene + 2)] = concNxp #Km
        rxnWs['I%d' % (iGene + 2)] = 'EC' #Cross reference source
        rxnWs['J%d' % (iGene + 2)] = '2.7.7.6' #Cross reference ID
        
        #translation reaction
        aaCnts = []
        for aa in aminoAcids:
            aaCnts.append('(%d) %s' % (protSeqs[iGene].seq.count(aa['id']), aa['name'], ))
            
        rxnWs['A%d' % (iGene + 2 + nGene)] = 'Translation-%s' % gene['id'] #ID
        rxnWs['B%d' % (iGene + 2 + nGene)] = 'Translation (%s)' % gene['name'] #name
        rxnWs['C%d' % (iGene + 2 + nGene)] = 'Translation' #submodel
        rxnWs['D%d' % (iGene + 2 + nGene)] = '[c]: %s + (%d) GTP + (%d) H2O ==> %s-Protein + (%d) GDP + (%d) PI + (%d) H' % ( #stoichometry
            ' + '.join(aaCnts),
            2 * len(protSeqs[iGene].seq) + 3,
            2 * len(protSeqs[iGene].seq) + 3 - (len(protSeqs[iGene].seq) - 1),
            gene['id'],
            2 * len(protSeqs[iGene].seq),
            2 * len(protSeqs[iGene].seq),
            2 * len(protSeqs[iGene].seq),)
        rxnWs['E%d' % (iGene + 2 + nGene)] = 'Ribosome-Protein' #enzyme
        rxnWs['F%d' % (iGene + 2 + nGene)] = 'Vmax * min(Ala[c], Arg[c], Asp[c], Asn[c], Cys[c], Gln[c], Glu[c], Gly[c], His[c], Ile[c], Leu[c], Lys[c], Met[c], Phe[c], Pro[c], Ser[c], Thr[c], Trp[c], Tyr[c], Val[c]) / (Km + min(Ala[c], Arg[c], Asp[c], Asn[c], Cys[c], Gln[c], Glu[c], Gly[c], His[c], Ile[c], Leu[c], Lys[c], Met[c], Phe[c], Pro[c], Ser[c], Thr[c], Trp[c], Tyr[c], Val[c])) * %s-Rna[c] * Ribosome-Protein[c]' % gene['id'] #rate law
        rxnWs['G%d' % (iGene + 2 + nGene)] = math.log(2) / CellCycleLength * 2 / nRnaCopy #Vmax
        rxnWs['H%d' % (iGene + 2 + nGene)] = concAa #Km
        rxnWs['I%d' % (iGene + 2 + nGene)] = 'EC' #Cross reference source
        rxnWs['J%d' % (iGene + 2 + nGene)] = '2.3.2.12' #Cross reference ID
        
        #RNA degradation reaction
        rxnWs['A%d' % (iGene + 2 + nGene * 2)] = 'RnaDegradation-%s' % gene['id'] #ID
        rxnWs['B%d' % (iGene + 2 + nGene * 2)] = 'RNA degradation (%s)' % gene['name'] #name
        rxnWs['C%d' % (iGene + 2 + nGene * 2)] = 'RnaDegradation' #submodel
        rxnWs['D%d' % (iGene + 2 + nGene * 2)] = '[c]: %s-Rna + (%d) H2O ==> (%d) AMP + (%d) CMP + (%d) GMP + (%d) UMP + (%d) H' % ( #stoichometry
            gene['id'],
            len(rnaSeqs[iGene].seq) - 1,
            rnaSeqs[iGene].seq.count('A'), 
            rnaSeqs[iGene].seq.count('C'), 
            rnaSeqs[iGene].seq.count('G'),
            rnaSeqs[iGene].seq.count('U'),
            len(rnaSeqs[iGene].seq) - 1, )
        rxnWs['E%d' % (iGene + 2 + nGene * 2)] = 'Rnase-Protein' #enzyme
        rxnWs['F%d' % (iGene + 2 + nGene * 2)] = 'Vmax * %s-Rna[c] / (%s-Rna[c] + Km) * Rnase-Protein[c]' % (gene['id'], gene['id']) #rate law
        rxnWs['G%d' % (iGene + 2 + nGene * 2)] = math.log(2) / RnaHalfLife * 2 * nRnaCopy / nProtCopy #Vmax
        rxnWs['H%d' % (iGene + 2 + nGene * 2)] = nRnaCopy / cellVol / nAvogadro #Km
        rxnWs['I%d' % (iGene + 2 + nGene * 2)] = 'EC' #Cross reference source
        rxnWs['J%d' % (iGene + 2 + nGene * 2)] = '3.1.-.-' #Cross reference ID
        
    #biomass production
    baseCntsNtp = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    nNtp_trn = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    nNmp_dcy = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    nAa_trl = {}
    for aa in aminoAcids:
        nAa_trl[aa['id']] = 0
    for iGene, gene in enumerate(genes):
        baseCntsNtp['A'] += rnaSeqs[iGene].seq.count('A') * nRnaCopy
        baseCntsNtp['C'] += rnaSeqs[iGene].seq.count('C') * nRnaCopy
        baseCntsNtp['G'] += rnaSeqs[iGene].seq.count('G') * nRnaCopy
        baseCntsNtp['U'] += rnaSeqs[iGene].seq.count('U') * nRnaCopy
        for aa in aminoAcids:
            nAa_trl[aa['id']] += protSeqs[iGene].seq.count(aa['id']) * nProtCopy
                
    nNtp_trn['A'] = baseCntsNtp['A'] * (1 + CellCycleLength / RnaHalfLife)
    nNtp_trn['C'] = baseCntsNtp['C'] * (1 + CellCycleLength / RnaHalfLife)
    nNtp_trn['G'] = baseCntsNtp['G'] * (1 + CellCycleLength / RnaHalfLife)
    nNtp_trn['U'] = baseCntsNtp['U'] * (1 + CellCycleLength / RnaHalfLife)
    
    nNmp_dcy['A'] = baseCntsNtp['A'] * CellCycleLength / RnaHalfLife
    nNmp_dcy['C'] = baseCntsNtp['C'] * CellCycleLength / RnaHalfLife
    nNmp_dcy['G'] = baseCntsNtp['G'] * CellCycleLength / RnaHalfLife
    nNmp_dcy['U'] = baseCntsNtp['U'] * CellCycleLength / RnaHalfLife
    
    nPPi_trn = (protLen + 1) * 3 * nRnaCopy * len(genes) * (1 + CellCycleLength / RnaHalfLife)
    nH2O_trn = nRnaCopy * len(genes) * (1 + CellCycleLength / RnaHalfLife)
    nH2O_dcy = ((protLen + 1) * 3 - 1) * nRnaCopy * len(genes) * (CellCycleLength / RnaHalfLife)
   
    nGtp_trl = (2 * protLen + 3) * nProtCopy * len(genes)
    nH2o_trl = (protLen + 2) * nProtCopy * len(genes)
        
    stochArr = []
    stochArr.append({'species': 'Biomass', 'coeff': 1})
    
    stochArr.append({'species': 'H2O', 'coeff': 0
        - cellVol * concH2O * nAvogadro 
        - nH2O_trn 
        - nH2o_trl
        - nH2O_dcy
        })
    stochArr.append({'species': 'H', 'coeff': 0
        - cellVol * concH * nAvogadro
        + nH2O_trn
        + nGtp_trl
        + nH2O_dcy
        })
        
    stochArr.append({'species': 'ATP', 'coeff': -cellVol * concNxp * nAvogadro - nNtp_trn['A']})
    stochArr.append({'species': 'CTP', 'coeff': -cellVol * concNxp * nAvogadro - nNtp_trn['C']})
    stochArr.append({'species': 'GTP', 'coeff': -cellVol * concNxp * nAvogadro - nNtp_trn['G'] - nGtp_trl})
    stochArr.append({'species': 'UTP', 'coeff': -cellVol * concNxp * nAvogadro - nNtp_trn['U']})
    
    stochArr.append({'species': 'GDP', 'coeff': -cellVol * concNxp * nAvogadro + nGtp_trl})
    
    stochArr.append({'species': 'AMP', 'coeff': -cellVol * concNxp * nAvogadro + nNmp_dcy['A']})
    stochArr.append({'species': 'CMP', 'coeff': -cellVol * concNxp * nAvogadro + nNmp_dcy['C']})
    stochArr.append({'species': 'GMP', 'coeff': -cellVol * concNxp * nAvogadro + nNmp_dcy['G']})
    stochArr.append({'species': 'UMP', 'coeff': -cellVol * concNxp * nAvogadro + nNmp_dcy['U']})
    
    stochArr.append({'species': 'PPI', 'coeff': -cellVol * concPPi * nAvogadro + nPPi_trn})
    stochArr.append({'species': 'PI', 'coeff': -cellVol * concPi * nAvogadro + nGtp_trl})
    
    for aa in aminoAcids:
        stochArr.append({'species': aa['name'], 'coeff': -cellVol * concAa * nAvogadro - nAa_trl[aa['id']]})
    
    lhs = []
    rhs = []
    for part in stochArr:
        if part['coeff'] < 0:
            lhs.append('(%e) %s' % (-part['coeff'], part['species']))
        else:
            rhs.append('(%e) %s' % ( part['coeff'], part['species']))
    
    stoichiometry = '[c]: %s ==> %s' % (' + '.join(lhs), ' + '.join(rhs))
        
    rxnWs['A%d' % (nGene * 3 + 2)] = 'BiomassProduction' #id
    rxnWs['B%d' % (nGene * 3 + 2)] = 'Biomass production' #id
    rxnWs['C%d' % (nGene * 3 + 2)] = 'Metabolism' #id
    rxnWs['D%d' % (nGene * 3 + 2)] = stoichiometry #stoichiometry
    rxnWs['E%d' % (nGene * 3 + 2)] = '' #enzyme
    rxnWs['F%d' % (nGene * 3 + 2)] = '' #rate law
    rxnWs['G%d' % (nGene * 3 + 2)] = '' #Vmax
    rxnWs['H%d' % (nGene * 3 + 2)] = '' #Km
    rxnWs['I%d' % (nGene * 3 + 2)] = '' #Cross reference source
    rxnWs['J%d' % (nGene * 3 + 2)] = '' #Cross reference id
    
    #style
    specWs.freeze_panes = specWs['B2']
    rxnWs.freeze_panes = rxnWs['B2']   
    
    for i in range(1, 12 + 1):
        cell = specWs.cell(row=1, column=i)
        cell.font = cell.font.copy(bold = True)
        cell.fill = cell.fill.copy(fgColor=Color('CCCCCC'), patternType='solid')
        
    for i in range(1, 10 + 1):
        cell = rxnWs.cell(row=1, column=i)
        cell.font = cell.font.copy(bold = True)
        cell.fill = cell.fill.copy(fgColor=Color('CCCCCC'), patternType='solid')
        
    for i in range(1, 2 * nGene + 2):
        specWs.row_dimensions[i].height = 15
        for j in range(1, 12 + 1):
            cell = specWs.cell(row=i, column=j)
            cell.alignment = cell.alignment.copy(wrap_text = True)
            
    for i in range(1, 3 * nGene + 3):
        rxnWs.row_dimensions[i].height = 15
        for j in range(1, 10 + 1):
            cell = rxnWs.cell(row=i, column=j)
            cell.alignment = cell.alignment.copy(wrap_text = True)
    
    wb.save("Model-RNA and proteins.xlsx")    

    '''save RNA and protein sequences in FASTA format'''
    output_handle = open("Genes.fasta", "w")
    SeqIO.write(rnaSeqs, output_handle, "fasta")
    output_handle.close()
    
    output_handle = open("RNAs.fasta", "w")
    SeqIO.write(rnaSeqs, output_handle, "fasta")
    output_handle.close()
    
    output_handle = open("Proteins.fasta", "w")
    SeqIO.write(protSeqs, output_handle, "fasta")
    output_handle.close()
    
#Generates random protein-coding DNA sequence of length 3 * (protLen + 1) that begins with a start codon and ends with a stop codon such that the corresponding protein has length protLen.
def generateRandomSequence(protLen = 100, startCodon = 'ATG', stopCodon = 'TAG'):
    dnaBases = 'ACGT'
    
    dnaSeqStr = startCodon
    for iCodon in range(protLen-1):
        validCodon = False
        while not(validCodon):
            codonStr = '' \
                + dnaBases[random.randint(0, 4)] \
                + dnaBases[random.randint(0, 4)] \
                + dnaBases[random.randint(0, 4)]
            validCodon = Seq(dnaSeqStr + codonStr).translate(table=4).find('*') == -1
        dnaSeqStr += codonStr
    dnaSeqStr += stopCodon
        
    return Seq(dnaSeqStr)