"""ESM-2 650M cross-validation of the SLC5A6 missense pathogenicity analysis.
Independent check that the 'passenger' conclusion is robust to model scale.
Reuses exact parsing/BLOSUM logic from esm_pathogenicity_mc3.py; only the model changes.
"""
import os
import re
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

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

AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
       "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
       "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}

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

def parse_missense(maf_path):
    maf = pd.read_csv(maf_path, sep='\t', dtype=str)
    ms = maf[maf['Variant_Classification'] == 'Missense_Mutation'].copy()
    muts = []
    for _, r in ms.iterrows():
        h = r['HGVSp_Short']
        if not isinstance(h, str) or h == '.':
            continue
        h = h.split('/')[0]
        m = re.match(r'p\.([A-Za-z]{3}|[A-Z])(\d+)([A-Za-z]{3}|[A-Z])', h.strip())
        if not m:
            continue
        r1, pos, a1 = m.group(1), int(m.group(2)), m.group(3)
        ref = AA3.get(r1, r1 if len(r1) == 1 else '?')
        alt = AA3.get(a1, a1 if len(a1) == 1 else '?')
        if ref == '?' or alt == '?' or pos < 1 or pos > 635:
            continue
        muts.append({'pos': pos, 'ref': ref, 'alt': alt, 'hgvsp': h})
    seen, uniq = set(), []
    for m in muts:
        k = (m['pos'], m['ref'], m['alt'])
        if k not in seen:
            seen.add(k); uniq.append(m)
    return uniq

def main():
    uniq = parse_missense('D:/Researching/SMVT/02_Data/cleaned/slc5a6_mc3_maf.csv')
    uniq.append({'pos': 123, 'ref': 'R', 'alt': 'L', 'hgvsp': 'p.R123L', 'source': 'ClinVar'})
    for m in uniq:
        m['domain'] = domain_of(m['pos'])
    print(f'unique missense + control: {len(uniq)}')

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model_name = 'facebook/esm2_t33_650M_UR50D'
    print(f'loading {model_name} on {device}...')
    tok = AutoTokenizer.from_pretrained(model_name, clean_up_tokenization_spaces=True)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    model = model.to(device); model.eval()

    results = []
    with torch.no_grad():
        for i, m in enumerate(uniq):
            seq = list(SMVT_SEQUENCE)
            seq[m['pos'] - 1] = '<mask>'
            inputs = tok(''.join(seq), return_tensors='pt').to(device)
            out = model(**inputs).logits[0]
            mask_idx = inputs['input_ids'][0].tolist().index(tok.mask_token_id)
            probs = torch.softmax(out[mask_idx], dim=-1).cpu().numpy()
            wt_idx = tok.convert_tokens_to_ids(m['ref'])
            mt_idx = tok.convert_tokens_to_ids(m['alt'])
            llr = float(np.log(probs[mt_idx]) - np.log(probs[wt_idx])) if wt_idx < len(probs) and mt_idx < len(probs) else np.nan
            results.append({**m, 'esm650_llr': llr})
            if i % 20 == 0:
                print(f'  {i}/{len(uniq)}')

    df = pd.DataFrame(results)
    df['esm650_percentile'] = df['esm650_llr'].rank(pct=True) * 100
    df = df.sort_values('esm650_llr')
    out = 'D:/Researching/SMVT/03_Analysis/outputs/esm650_mc3.csv'
    if os.path.exists(out):
        import datetime
        out = out.replace('.csv', f'_{datetime.datetime.now().strftime("%H%M%S")}.csv')
    df.to_csv(out, index=False)
    print(f'\nDONE: {len(df)} rows -> {out}')
    print('Most destabilizing (650M):')
    print(df.head(10)[['hgvsp','pos','domain','esm650_llr','esm650_percentile']].to_string(index=False))
    ctl = df[df['hgvsp'] == 'p.R123L']
    if len(ctl):
        print(f'\nR123L control: llr={ctl.iloc[0]["esm650_llr"]:.3f}, percentile={ctl.iloc[0]["esm650_percentile"]:.0f}')

if __name__ == '__main__':
    main()
