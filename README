# LoRa synchronization test source code
TODO: paper reference

## Description
*run_test.py* tool generates a .csv file containing results of the proposed synchronization 
algorithm and the exhaustive algorithm on randomly generated LoRa frames.

Program arguments:
```
    --sf' : chossen spreading factor, default 7
    --snr : signal noise ratio, default 0
    --count': number of frame tested, default 1000
    --output-file: output file name, default='results.csv'
```

Run:
```
    $ python3 run_test.py --sf 7 --snr 0 --count 1000 --output-file results.csv
```

The output csv contains the following columns:
- Date : date of the test
- SF : spreading factor
- SNR : signal noise ratio
- Delay : random signal delay
- Offset : random signal frequency offset
- ExDelay : exhaustive synchronization estimated delay
- ExOffset : exhaustive synchronization estimated offset
- SynDelay : proposed synchronization estimated delay
- SynOffset : proposed synchronization estimated offset

## Install prerequisite

### Manual package installation

python: version 3.5 or newer
required package : numpy

### Using pipenv

- Install python 3.7
- Install pipenv with pip:
    $ pip3 install --user pipenv

- Install dependencies
    $ pipenv --python 3.7
    $ pipenv shell
    $ pipenv install
