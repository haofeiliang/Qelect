



# Qelect: Lattice-based Single Secret Leader Election Made Practical

## Abstract:
In a single secret leader election (SSLE) protocol, all parties collectively and obliviously elect one leader. No one else should learn its identity unless it reveals itself as the leader.
The problem is first formalized by Boneh et al. (AFT'20), which proposes an efficient construction based on the Decision Diffie-Hellman (DDH) assumption.
Considering the potential risk of quantum computers, several follow-ups focus on designing a post-quantum secure SSLE protocol based on pure lattices or fully homomorphic encryption. However, no concrete benchmarks demonstrate the feasibility of deploying such heavy cryptographic primitives.

In this work, we present Qelect, the first practical constant-round post-quantum secure SSLE protocol.
We first adapt the commitment scheme in Boneh et al. (AFT'23) into a multi-party randomizable commitment scheme, and propose our novel construction based on an adapted version of ring learning with errors (RLWE) problem.
We then use it as a building block and construct a constant-round single secret leader election (crSSLE) scheme.
We utilize the single instruction multiple data (SIMD) property of a specific threshold fully homomorphic encryption (tFHE) scheme to evaluate our election circuit efficiently.
Finally, we built Qelect from the crSSLE scheme, with performance optimizations including a preprocessing phase to amortize the local computation runtime and a retroactive detection phase to avoid the heavy zero-knowledge proofs during the election phase.
Qelect achieves asymptotic improvements and is concretely practical.
We implemented a prototype of Qelect and evaluated its performance in a WAN.
Qelect is at least two orders of magnitude faster than the state-of-the-art.

## License

The Qelect library is developed by [Yunhao Wang](https://scholar.google.com/citations?user=-3s-pjIAAAAJ&hl=en) and [Fan Zhang](https://www.fanzhang.me/), and is released under the MIT License (see the LICENSE file).

## Descriptions of library
The main simulation code is contained in ```MPE_docker_version.cpp```, with ```regevEncryption.h``` including the assistant functions for (R)LWE encryption scheme  and ```util.h``` including all other assistant functions.
```global.h``` is used to share common parameters among scripts.
```execute_commands.sh``` is used to simulate the communication between multiple instances when analyzing the communication time under LAN/WAN setting.


## The instructions to run the project via Dockerfile:

### install docker
```
sudo apt-get update
sudo apt-get install -y ca-certificates curl
# might not needed, ignore for now: sudo apt-get install -y -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$UBUNTU_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

### execute the qelect
```
# download the zip file and unzip the folder under MPS
git clone --depth 1 https://github.com/haofeiliang/Qelect.git qelect
cd qelect

# build without cache
sudo docker build --network host --no-cache -t qelect_project .

# build with proxy
sudo docker build --network host \
      --build-arg HTTP_PROXY=http://localhost:10808 \
      --build-arg HTTPS_PROXY=http://localhost:10808 \
      -t qelect_project .

# run
sudo docker run qelect_project
```

An example output is as follows:
```
After param generation.
Initial noise budget: 613
Lefted noise budget: 47
Total broadcast communication size (for each party): 2147 MB.
---------------------------------
Total number of parties: 1024
Preprocessed time      : 136201928 us.
Total time             : 5921790 us.
```
The local computation time for a single party, with party size = 1024, is ~5.9 seconds.

## The instructions to run the project via cmake:
```
# If permission required, please add sudo before the commands as needed

sudo apt-get update && sudo apt-get install build-essential # if needed
sudo apt-get install autoconf # if no autoconf
sudo apt-get install cmake # if no cmake
sudo apt-get install libgmp3-dev # if no gmp
sudo apt-get install libntl-dev=11.4.3-1build1 # if no ntl, or sudo apt-get install libntl-dev is version not found

# change build_path to where you want the dependency libraries installed
MPSDIR=~
BUILDDIR=$MPSDIR/MPS/build

# used for data type definitions and distribution, generators, etc
cd $MPSDIR && git clone -b v1.11.9 https://gitlab.com/palisade/palisade-release
cd palisade-release
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$BUILDDIR
make
make install


# a separate fork of SEAL library that overwrite some private functions, used in prior works
# we also depend on this library which is used by the prior work we are based on
cd $MPSDIR && git clone https://github.com/wyunhao/SEAL
cd SEAL
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$BUILDDIR -DSEAL_USE_INTEL_HEXL=ON 
cmake --build build
cmake --install build


cd $MPSDIR/MPS/build
cmake .. -DCMAKE_PREFIX_PATH=$BUILDDIR
mkdir ../data
mkdir ../data/perm
make
```
