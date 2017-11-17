from benchmark import StandardRunner

for version in ['multi','one']:
    sr = StandardRunner()
    sr.run_standard_flex(version=version, name_prefix='flex'+version+'chain_')