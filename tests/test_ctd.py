import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class CTDPreparation(unittest.TestCase):
    def test_species_deduplication_and_sampling_order(self):
        spec = importlib.util.spec_from_file_location('prepare_ctd', ROOT / 'tools/prepare_ctd.py')
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        def row(chemical, gene, species):
            return '\t'.join([chemical,'ID','CAS',gene,'GID','','organism',species,'effect','action','reference'])
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'source.tsv.gz'
            target = Path(tmp) / 'pairs.tsv'
            with gzip.open(source,'wt',encoding='utf-8') as f:
                f.write('# synthetic CTD fixture\n')
                f.write('\n'.join([row('Zeta','GENE2','9606'),row('Alpha','GENE1','9606'),
                                   row('Alpha','GENE1','9606'),row('Other','GENE3','10090')])+'\n')
            self.assertEqual(mod.extract(source,target),2)
            self.assertEqual(target.read_text(encoding='utf-8'),'Alpha\tGENE1\nZeta\tGENE2\n')

if __name__ == '__main__':
    unittest.main()
