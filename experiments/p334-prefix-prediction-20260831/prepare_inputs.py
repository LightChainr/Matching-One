#!/usr/bin/env python3
"""Extract sufficient paired prefix statistics; do not fit a prediction model."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

SOURCES={
    "results-extension/prefix_statistics_N325.npz":"be9d89e8f2eb7bf8cb9ce2168ac8737dc1466ac89d6540f8250ebf98014ef79e",
    "results-extension/prefix_statistics_N425.npz":"fa79cf810c01bca1cf560f0d1c53cf0bd47e665f9acbb1a33dd45a16f39c3154",
    "results-exact-score/prefix_statistics_N325.npz":"0a8bd345c603a3583075b84649f5c0c2454af2f055275cb4574dea1a7c085492",
    "results-exact-score/prefix_statistics_N425.npz":"65d298c3331ab77a26d34ae19df188608a919074bca5bb462cb6488ef3efe957",
}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path(__file__).parent.parent/'p334-mechanism-response-20260831')
    p.add_argument('--output',type=Path,default=Path(__file__).parent/'inputs')
    a=p.parse_args();a.output.mkdir(exist_ok=False)
    manifest={'repository':'https://github.com/LightChainr/Matching-One','source_commit':'8ad30617b0a3076a5c01a208eb213096d8879b32',
              'source_package':'experiments/p334-mechanism-response-20260831','sources':[], 'derived_files':[], 'checks':{}}
    for path,digest in SOURCES.items():
        assert hashlib.sha256((a.source/path).read_bytes()).hexdigest()==digest,path
        manifest['sources'].append({'path':path,'sha256':digest})
    for n in (325,425):
        full=np.load(a.source/f'results-exact-score/prefix_statistics_N{n}.npz')
        extra=np.load(a.source/f'results-extension/prefix_statistics_N{n}.npz')
        choose=full['cell']==0
        assert np.array_equal(full['counter'][choose],extra['counter'])
        assert np.array_equal(full['batch'][choose],extra['batch'])
        fi={str(k):i for i,k in enumerate(full['labels'])}
        ei={str(k):i for i,k in enumerate(extra['labels'])}
        gram=np.stack([full['values'][choose,fi[f'exact_score_G.matrix[{o},{s}]']] for o in ('first','second') for s in ('first','second')],axis=-1).reshape(-1,2,2)
        assert np.max(abs(gram-gram.transpose(0,2,1)))<1e-15
        assert np.all(np.linalg.eigvalsh(gram)>0)
        means,energies=[],[]
        mean_error=energy_error=0.
        for obs,scale in (('p_ref.A',1.),('p_integral.A',-.5)):
            ix=[ei[f'{obs}.mean_J[{o},{s}]'] for o in ('first','second') for s in ('first','second')]
            j8,j64,j72=(extra[k][:,ix].reshape(-1,2,2) for k in ('old8','new64','combined72'))
            ee=ei[f'{obs}.E_frobenius_JZ_squared']
            u8,u64,u72=(extra[k][:,ee] for k in ('old8','new64','combined72'))
            reconstructed_mean=(8*j8+64*j64)/72
            reconstructed_energy=(8*7*u8+64*63*u64+2*8*64*np.sum(j8*j64,axis=(1,2)))/(72*71)
            mean_error=max(mean_error,float(np.max(abs(j72-reconstructed_mean))))
            energy_error=max(energy_error,float(np.max(abs(u72-reconstructed_energy))))
            assert np.allclose(j72,reconstructed_mean,rtol=1e-11,atol=1e-16)
            assert np.allclose(u72,reconstructed_energy,rtol=1e-10,atol=1e-18)
            means.append(scale*j72);energies.append(scale*scale*u72)
        target=a.output/f'N{n}.npz'
        np.savez_compressed(target,N=n,counter=extra['counter'],batch=extra['batch'],fold=extra['batch']%5,
                            Jbar=np.stack(means,axis=1),energy_U=np.stack(energies,axis=1),G=gram,
                            output_names=np.array(['A_p_ref','C_over_Nplus1']),quartets=72,
                            cell=np.zeros(len(gram),dtype=np.int64))
        manifest['derived_files'].append({'path':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'prefixes':len(gram)})
        manifest['checks'][str(n)]={'cell00_prefixes':len(gram),'matrix_mean_combination_max_error':mean_error,'energy_combination_max_error':energy_error,
                                    'gram_eigenvalue_min':float(np.linalg.eigvalsh(gram).min()),'original_batch_counts':np.bincount(extra['batch'],minlength=20).tolist()}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest['checks'],indent=2))

if __name__=='__main__':main()
