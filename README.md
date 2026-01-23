![Elle logo](elle.jpg) 
# elle-linux
linux compatibly version of elle   


This repository is the **release, modern Linux version** of the ELLE software, adapted and maintained by **Maxime Fatzaun**.  
It is a curated mirror of the private development repository and may be updated periodically.

> **Status — 2025-09-29:** Adapted and verified to compile on modern Linux distributions.

## What this is
- Linux-focused, modernized build of ELLE.
- Public release mirror; issues and PRs are welcome here. Accepted changes may be tested in the private dev repo before merging.

## Contact
Questions or queries? Contact **Maxime Fatzaun** at **<maxime.fatzaun@fau.de>**.

## License
**ELLE Project — License Notice**

This repository includes legacy ELLE geological simulation sources and related materials.  
The primary license text is provided in:

elle/elle/COPYING.txt   

Please refer to that file for the full terms and conditions. If multiple components
carry their own licenses, they are located alongside their sources within the `elle/` tree.   

Elle Microstructure Modelling Package Version 2     
Copyright (C) 2006  Lynn Evans, Jens Becker, Mark Jessell, Daniel Koehn et al   

### Notes from Daniel
- To install Elle, execute the script **install_elle.sh** in the **scripts directory**, then everything will be installed and compiled. To make the script executable, if necessary, use the command: chmod +x install_elle.sh. 

- The install script is for Linux, and I tested it on Ubuntu 22.02. There are major differences between Linux versions, so older versions require adjustment, as the one on the desktop computer in my office.

- If there are any issues, bugs, problems, or difficulty installing the software or if you have questions, feel free to contact Daniel Koehn.

Also, to recompile all of ELLE after the install, to take the code modifications into account, use: cmake —build build -j

I installed it on Ubuntu 24 as well 
It worked without any issue on one computer but on the second one I had to install this by hand: 

wget https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
tar -zxvf gsl-2.7.tar.gz
cd gsl-2.7
./configure
make
sudo make install
sudo ldconfig
