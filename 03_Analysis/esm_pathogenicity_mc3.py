"""ESM-2 pathogenicity for 96 SLC5A6 MC3 missense mutations + germline controls.
Fixes the stale 'zero missense' analysis that used incomplete cBioPortal data."""
import os
import pandas as pd
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

# ── SMVT Sequence (UniProt Q9Y289) ──
SMVT_SEQUENCE = (
    "MSVGVSTSAPLSPTSGTSVGMSTFSIMDYVVFVLLLVLSLAIGLYHACRGWGRHTVGELL"
    "MADRKMGCLPVALSLLATFQSAVAILGVPSEIYRFGTQYWFLGCCYFLGLLIPAHIFIPV"
    "FYRLHLTSAYEYLELRFNKTVRVCGTVTFIFQMVIYMGVVLYAPSLALNAVTGFDLWLSV"
    "LALGIVCTVYTALGGLKAVIWTDVFQTLVMFLGQLAVIIVGSAKVGGLGRVWAVASQHGR"
    "ISGFELDPDPFVRHTFWTLAFGGVFMMLSLYGVNQAQVQRYLSSRTEKAAVLSCYAVFPF"
    "QQVSLCVGCLIGLVMFAYYQEYPMSIQQAQAAPDQFVLYFVMDLLKGLPGLPGLFIACLF"
    "SGSLSTISSAFNSLATVTMEDLIRPWFPEFSEARAIMLSRGLAFGYGLLCLGMAYISSQM"
    "GPVLQAAISIFGMVGGPLLGLFCLGMFFPCANPPGAVVGLLAGLVMAFWIGIGSIVTSMG"
    "SSMPPSPSNGSSFSLPTNLTVATVTTLMPLTTFSKPTGLQRFYSLSYLWYSAHNSTTVIV"
    "VGLIVSLLTGRMRGRSLNPATIYPVLPKLLSLLPLSCQKRLHCRSYGQDHLDTGLFPEKP"
    "RNGVLGDSRDKEAMALDGTAYQGSSSTCILQETSL"
)
assert len(SMVT_SEQUENCE) == 635

DOMAINS = {
    "N-term": (1, 80),
    "Transmembrane_helix_1": (81, 115),
    "Transmembrane_helix_2": (120, 155),
    "Transmembrane_helix_3": (165, 195),
    "Transmembrane_helix_4": (200, 235),
    "Transmembrane_helix_5": (245, 275),
    "Transmembrane_helix_6": (280, 315),
    "Transmembrane_helix_7": (325, 360),
}
def domain_of(pos):
    for name, (s, e) in DOMAINS.items():
        if s <= pos <= e:
            return name
    return "Loop"

AA1 = "ACDEFGHIKLMNPQRSTVWY"
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
       "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
       "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
AA1_TO_3 = {v:k for k,v in AA3.items()}

# BLOSUM62 (reduced to 20 AA)
BLOSUM62 = {
 ('A','A'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,('A','Q'):-1,('A','E'):-1,('A','G'):0,('A','H'):-2,('A','I'):-1,('A','L'):-1,('A','K'):-1,('A','M'):-1,('A','F'):-2,('A','P'):-1,('A','S'):1,('A','T'):0,('A','W'):-3,('A','Y'):-2,('A','V'):0,
 ('R','R'):5,('R','N'):0,('R','D'):-2,('R','C'):-3,('R','Q'):1,('R','E'):0,('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,('R','K'):2,('R','M'):-1,('R','F'):-3,('R','P'):-2,('R','S'):-1,('R','T'):-1,('R','W'):-3,('R','Y'):-2,('R','V'):-3,
 ('N','N'):6,('N','D'):1,('N','C'):-3,('N','Q'):0,('N','E'):0,('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,('N','K'):0,('N','M'):-2,('N','F'):-3,('N','P'):-2,('N','S'):1,('N','T'):0,('N','W'):-4,('N','Y'):-2,('N','V'):-3,
 ('D','D'):6,('D','C'):-3,('D','Q'):0,('D','E'):2,('D','G'):-1,('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,('D','M'):-3,('D','F'):-3,('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,('D','Y'):-3,('D','V'):-3,
 ('C','C'):9,('C','Q'):-3,('C','E'):-4,('C','G'):-3,('C','H'):-3,('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,
 ('Q','Q'):5,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,('Q','L'):-2,('Q','K'):1,('Q','M'):0,('Q','F'):-3,('Q','P'):-1,('Q','S'):0,('Q','T'):-1,('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,
 ('E','E'):5,('E','G'):-2,('E','H'):0,('E','I'):-3,('E','L'):-3,('E','K'):1,('E','M'):-2,('E','F'):-3,('E','P'):-1,('E','S'):0,('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,
 ('G','G'):6,('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,('G','M'):-3,('G','F'):-3,('G','P'):-2,('G','S'):0,('G','T'):-2,('G','W'):-2,('G','Y'):-3,('G','V'):-3,
 ('H','H'):8,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,('H','F'):-1,('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,('H','Y'):2,('H','V'):-3,
 ('I','I'):4,('I','L'):2,('I','K'):-3,('I','M'):1,('I','F'):0,('I','P'):-3,('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,
 ('L','L'):4,('L','K'):-2,('L','M'):2,('L','F'):0,('L','P'):-3,('L','S'):-2,('L','T'):-1,('L','W'):-2,('L','Y'):-1,('L','V'):1,
 ('K','K'):5,('K','M'):-1,('K','F'):-3,('K','P'):-1,('K','S'):0,('K','T'):-1,('K','W'):-3,('K','Y'):-2,('K','V'):-2,
 ('M','M'):5,('M','F'):-1,('M','P'):-2,('M','S'):-1,('M','T'):-1,('M','W'):-1,('M','Y'):-1,('M','V'):1,
 ('F','F'):6,('F','P'):-4,('F','S'):-2,('F','T'):-2,('F','W'):1,('F','Y'):3,('F','V'):-1,
 ('P','P'):7,('P','S'):-1,('P','T'):-2,('P','W'):-4,('P','Y'):-3,('P','V'):-2,
 ('S','S'):4,('S','T'):1,('S','W'):-3,('S','Y'):-2,('S','V'):-2,
 ('T','T'):5,('T','W'):-2,('T','Y'):-2,('T','V'):0,
 ('W','W'):11,('W','Y'):2,('W','V'):-3,
 ('Y','Y'):7,('Y','V'):-1,
 ('V','V'):4,
}
def blosum62(a, b):
    return BLOSUM62.get(tuple(sorted((a, b))), BLOSUM62.get((a, b), -1))

def main():
    maf = pd.read_csv('D:/Researching/SMVT/02_Data/cleaned/slc5a6_mc3_maf.csv', sep='\t', dtype=str)
    ms = maf[maf['Variant_Classification'] == 'Missense_Mutation'].copy()
    print(f'MC3 missense: {len(ms)}')

    # parse HGVSp_Short -> (pos, ref, alt)
    muts = []
    for _, r in ms.iterrows():
        h = r['HGVSp_Short']
        if not isinstance(h, str) or h == '.':
            continue
        # formats: p.L461M / p.Val120Ile / p.V120I
        h = h.split('/')[0]
        import re
        m = re.match(r'p\.([A-Za-z]{3}|[A-Z])(\d+)([A-Za-z]{3}|[A-Z])', h.strip())
        if not m:
            continue
        r1, pos, a1 = m.group(1), int(m.group(2)), m.group(3)
        ref = AA3.get(r1, r1 if len(r1) == 1 else '?')
        alt = AA3.get(a1, a1 if len(a1) == 1 else '?')
        if ref == '?' or alt == '?' or pos < 1 or pos > 635:
            continue
        muts.append({'pos': pos, 'ref': ref, 'alt': alt, 'hgvsp': h, 'patient': r['Tumor_Sample_Barcode'][:12], 'cancer': r['Cancer'] if 'Cancer' in maf.columns else ''})
    # dedupe by (pos, ref, alt)
    seen = set(); uniq = []
    for m in muts:
        k = (m['pos'], m['ref'], m['alt'])
        if k not in seen:
            seen.add(k); uniq.append(m)
    print(f'unique missense: {len(uniq)}')

    # germline controls
    germline = [{'pos': 123, 'ref': 'R', 'alt': 'L', 'source': 'ClinVar', 'note': 'biotin-responsive neurodegeneration'}]
    for g in germline:
        uniq.append({**g, 'patient': 'Germline', 'cancer': 'Germline'})

    # ESM-1v
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'loading ESM-2 on {device}...')
    tok = AutoTokenizer.from_pretrained('facebook/esm2_t12_35M_UR50D', revision='main', clean_up_tokenization_spaces=True)
    model = AutoModelForMaskedLM.from_pretrained('facebook/esm2_t12_35M_UR50D')
    model = model.to(device); model.eval()

    results = []
    with torch.no_grad():
        for i, m in enumerate(uniq):
            seq = list(SMVT_SEQUENCE)
            seq[m['pos'] - 1] = '<mask>'
            masked = ''.join(seq)
            inputs = tok(masked, return_tensors='pt').to(device)
            out = model(**inputs).logits[0]
            mask_idx = inputs['input_ids'][0].tolist().index(tok.mask_token_id)
            probs = torch.softmax(out[mask_idx], dim=-1).cpu().numpy()
            aa_idx = {tok.convert_ids_to_tokens(i): i for i in range(len(probs))}
            # LLR = log P(mutant) - log P(wildtype)
            wt_idx = tok.convert_tokens_to_ids(m['ref'])
            mt_idx = tok.convert_tokens_to_ids(m['alt'])
            llr = float(np.log(probs[mt_idx]) - np.log(probs[wt_idx])) if wt_idx < len(probs) and mt_idx < len(probs) else np.nan
            b62 = blosum62(m['ref'], m['alt'])
            results.append({**m, 'esm_llr': llr, 'blosum62': b62, 'domain': domain_of(m['pos'])})
            if i % 20 == 0:
                print(f'  {i}/{len(uniq)}')

    df = pd.DataFrame(results)
    df['esm_zscore'] = (df['esm_llr'] - df['esm_llr'].mean()) / df['esm_llr'].std()
    df['blosum_zscore'] = (df['blosum62'] - df['blosum62'].mean()) / df['blosum62'].std()
    df['composite'] = 0.7 * df['esm_zscore'] + 0.3 * df['blosum_zscore']
    df['composite_percentile'] = df['composite'].rank(pct=True) * 100
    df = df.sort_values('composite')
    out = 'D:/Researching/SMVT/03_Analysis/outputs/esm_pathogenicity_mc3.csv'
    if os.path.exists(out):
        import datetime
        out = out.replace('.csv', f'_{datetime.datetime.now().strftime("%H%M%S")}.csv')
    df.to_csv(out, index=False)
    print(f'\nDONE: {len(df)} rows -> {out}')
    print('\nMost pathogenic (lowest composite):')
    print(df.head(12)[['hgvsp','pos','ref','alt','domain','esm_llr','blosum62','composite','composite_percentile']].to_string())

if __name__ == '__main__':
    main()
