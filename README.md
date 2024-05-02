# UKBB

After cloning repo on UKBB virtual machine

1. Get cram list (using folder num)

```
./get_cram_list.sh ${fold_num}
```

2. install required programs

```
./unzip_make.sh
```

3. process crams as background jobs (done on per folder basis), will also upload finished files to UKKB

```
./parall_cram_run.sh ${fold_num}
```

5.
```
./tryptase.v1.R
```

