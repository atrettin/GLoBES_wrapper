# GLoBES_wrapper
Python wrapper for GLoBES
Author: Andrii Terliuk 


## Installation process

### Install GLoBES

First, load CVMFS software v. 4.0.1 into your environment!

Install GLoBES, currently available stable version is 3.0.11
Go to some folder and do: 
```
 wget http://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.0.11.tar.gz
 tar -xvzf globes-3.0.11.tar.gz
 cd ./globes-3.0.11
 ```
 
Install it using : 
```
./configure --prefix=<path/to/globes/installation/folder/>
make 
make install
```

### Install SNU
Get Sterile neutrinos and non-standard interactions extension for globes
https://www.mpi-hd.mpg.de/personalhomes/globes/tools.html

Place `snu.c` and `snu.h` into `./snu/` subfolder of this repository. 

go to `/snu` folder and apply a bit hacky patch to set correct electron density 
```
```
edit prefix path in Makefile:
```
prefix=<put/your/path/to/GLoBES/prefix/here>
```
and then do 
```
make
make install
```

### Install wrapper

After this you need to compile a Python interface for SNU extension. 
First, put location of your GLoBES/SNU extension installation as `prefix` to Makefile
Then run `make` and `make install`.
You should see `GLoBES.so` file and now it should be ready to use. 

## Testing the installation 

Go to examples folder and change path in `plot_numubar.py` file to the location of the folder containitng `GLoBES.so`. 
Then simply run the script, which should produce the following oscillogram: 


<object data="https://github.com/terliuk/GLoBES_wrapper/raw/master/examples/probability_example.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/terliuk/GLoBES_wrapper/raw/master/examples/probability_example.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/terliuk/GLoBES_wrapper/raw/master/examples/probability_example.pdf">Download PDF</a>.</p>
    </embed>
</object>
