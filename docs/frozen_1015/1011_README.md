# Publication analysis modules

`00_*` through `07_*` are the publication-facing analysis modules. All scripts
and third-party EasyColoc files in this tree are ordinary files or directories;
the tree contains no filesystem symlinks and can be transferred with `rsync`.

`SOURCE_MANIFEST.tsv` records where each materialized script set came from.
Module 02 `data/` and `outputs/` are intentionally empty contract directories:
controlled inputs and generated results were not copied into the code package.

The canonical entrances do not submit jobs automatically. Use the module-level
driver when one exists, and inspect its configuration before submitting.

Transfer the standalone publication tree without dereferencing options:

```bash
rsync -av --exclude='__pycache__/' 1011/ user@host:/path/to/apa-publication-code/
```

After transfer, run:

```bash
bash tests/check_portable_tree.sh
bash tests/smoke_module00.sh
bash tests/smoke_module06.sh
```
