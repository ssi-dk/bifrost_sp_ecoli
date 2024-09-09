# bifrost_sp_ecoli
This component runs given a sample_id already added into the bifrostDB. If the sample is registered as *Escherichia coli*, it should pull the paired reads, the contigs, and do the typing according to [ecoli_fbi](https://github.com/ssi-dk/ecoli_fbi)

## How to launch
```bash
git clone https://github.com/ssi-dk/bifrost_sp_ecoli.git
cd bifrost_sp_ecoli
git submodule init
git submodule update
bash install.sh -i LOCAL
conda activate bifrost_sp_ecoli_vx.x.x
export BIFROST_INSTALL_DIR='/your/path/'
BIFROST_DB_KEY="/your/key/here/" python -m bifrost_sp_ecoli -h
```
