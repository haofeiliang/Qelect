#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <thread>
#include "seal/util/polyarithsmallmod.h"
#include "seal/seal.h"
#include "math.h"

using namespace seal::util;
using namespace std;
using namespace seal;

vector<vector<uint64_t>> expand_ring_vector(vector<uint64_t>& a, int n, int q) {
    vector<vector<uint64_t>> res(n);

    for (int cnt = 0; cnt < n; cnt++) {
        res[cnt].resize(n);
        int ind = 0;
        for (int i = cnt; i >= 0 && ind < n; i--) {
            res[cnt][ind] = a[i];
            ind++;
        }

        for (int i = n-1; i > cnt && ind < n; i--) {
            res[cnt][ind] = q - a[i];
            ind++;
        }
    }

    return res;
}

vector<uint64_t> ring_multiply(vector<uint64_t>& a, vector<uint64_t>& b, int q, bool print = false) {
    int n = (int) a.size();
    vector<uint64_t> res(n);

    vector<vector<uint64_t>> expanded_A = expand_ring_vector(a, n, q);
    for (int i = 0; i < n; i++) {
        long temp = 0;
        for (int j = 0; j < n; j++) {
            temp = (temp + expanded_A[i][j] * b[j]) % q;
            temp = temp < 0 ? temp + q : temp;
        }
        res[i] = temp;
    }

    // if (print) {
    // for (int i = 0; i < n; i++) {
    //     cout << res[i] << ", " << res1[i] << "  ";
    // }
    // cout << endl;
    // }

    return res;
}


void saveCiphertext(Ciphertext& ct, const uint64_t index, const string folder = "perm/") {
    stringstream ss;
    ct.save(ss);
    ofstream datafile;
    datafile.open ("../data/"+folder+to_string(index)+".txt");
    datafile << ss.rdbuf();
    datafile.close();
}

void loadCiphertext(const SEALContext& context, Ciphertext& ct, const uint64_t index,
                    const string folder = "perm/") {

    // uint64_t ind = index / 2;
    ifstream datafile;
    datafile.open ("../data/"+folder+to_string(index)+".txt");

    stringstream buf;
    buf << datafile.rdbuf();

    ct.load(context, buf);
}

void saveDoubleVector(const vector<vector<uint64_t>>& input){
    ofstream datafile;
    datafile.open ("../data/perm/u_inverse.txt");

    for(size_t i = 0; i < input.size(); i++){
        for (size_t j = 0; j < input[0].size(); j++) {
            datafile << input[i][j] << "\n";
        }
    }

    datafile.close();
}

void loadDoubleVector(vector<vector<uint64_t>>& input){
    ifstream datafile;
    datafile.open ("../data/perm/u_inverse.txt");

    for(size_t i = 0; i < input.size(); i++){
        for (size_t j = 0; j < input[0].size(); j++) {
            uint64_t temp;
            datafile >> temp;
            input[i][j] = temp;
        }
    }

    datafile.close();
}


void print_ct_to_pl(Ciphertext& ct, SEALContext& context, SecretKey& sk, uint64_t ring_dim = 4) {
  Decryptor decryptor(context, sk);
  Plaintext pp;
  decryptor.decrypt(ct, pp);
  for (int i = 0; i < (int) ring_dim; i++) {
    cout << pp.data()[i] << " ";
  }
  cout << endl;
}


void print_ct_to_vec(Ciphertext& ct, SEALContext& context, SecretKey& sk, uint64_t ring_dim = 4) {
  Decryptor decryptor(context, sk);
  BatchEncoder batch_encoder(context);
  Plaintext pp;
  vector<uint64_t> msg(ring_dim);
  decryptor.decrypt(ct, pp);
  batch_encoder.decode(pp, msg);

  for (int i = 0; i < (int) ring_dim; i++) {
    cout << msg[i] << " ";
  }
  cout << endl;
}


Ciphertext extract_and_multiply_multi_core(const RelinKeys &relin_keys, const SEALContext& context,
                                           Ciphertext& ct1, Ciphertext& ct2, const int group_size = 4096,
                                           const int ring_dim = 32768, const int numcores = 8) {
    Evaluator evaluator(context);
    NTL::SetNumThreads(numcores);

    int sq_group_size = sqrt(group_size);
    int batch_size = ring_dim / sq_group_size; // extract each sq_group_size

    // vector<Ciphertext> add_tmp(sq_group_size);
    Ciphertext final_ct;
    Plaintext pl;
    pl.resize(ring_dim);
    pl.parms_id() = parms_id_zero;
    for (int j = 0; j < (int) ring_dim; j++) { 
        if (j < batch_size) {
            pl.data()[j] = 1;
        } else {
            pl.data()[j] = 0;
        }
    }
    Ciphertext tmp1, tmp2;
    evaluator.multiply_plain(ct1, pl, tmp1);
    evaluator.multiply_plain(ct2, pl, tmp2);
    evaluator.multiply(tmp1, tmp2, final_ct);

    int core_batch_size = batch_size / numcores;

    NTL_EXEC_RANGE(numcores, first, last);
    for (int c = first; c < last; c++) {
        for (int i = core_batch_size*c; i < core_batch_size*(c+1); i++) {

            if (i == 0) continue; // skip the initial setup

            // prepare the extraction plaintext
            Plaintext pl_core;
            pl_core.resize(ring_dim);
            pl_core.parms_id() = parms_id_zero;
            for (int j = 0; j < (int) ring_dim; j++) { 
                if (j >= sq_group_size*i && j < sq_group_size*(i+1)) {
                    pl_core.data()[j] = 1;
                } else {
                    pl_core.data()[j] = 0;
                }
            }

            Ciphertext tmp1_core, tmp2_core;
            evaluator.multiply_plain(ct1, pl_core, tmp1_core);
            evaluator.multiply_plain(ct2, pl_core, tmp2_core);

            Ciphertext tmp;
            evaluator.multiply(tmp1_core, tmp2_core, tmp);
            // evaluator.relinearize_inplace(tmp, relin_keys);
            evaluator.add_inplace(final_ct, tmp);
        }
    }
    NTL_EXEC_RANGE_END;

    evaluator.relinearize_inplace(final_ct, relin_keys);
    return final_ct;
}

Ciphertext extract_and_multiply(const RelinKeys &relin_keys, const SEALContext& context, Ciphertext& ct1,
              Ciphertext& ct2, const int group_size = 4096, const int ring_dim = 32768) {
    Evaluator evaluator(context);

    int sq_group_size = sqrt(group_size);
    int batch_size = ring_dim / sq_group_size; // extract each sq_group_size
    Plaintext pl;
    pl.resize(ring_dim);
    pl.parms_id() = parms_id_zero;
    // vector<Ciphertext> add_tmp(sq_group_size);
    Ciphertext final_ct;
    chrono::high_resolution_clock::time_point s1, e1;

    // s1 = chrono::high_resolution_clock::now();
    evaluator.transform_to_ntt_inplace(ct1);
    evaluator.transform_to_ntt_inplace(ct2);
    // e1 = chrono::high_resolution_clock::now();

    for (int i = 0; i < batch_size; i++) {
        // prepare the extraction plaintext
        for (int j = 0; j < (int) ring_dim; j++) { 
            // cout << "????? " << j << endl;
            if (j >= sq_group_size*i && j < sq_group_size*(i+1)) {
                pl.data()[j] = 1;
            } else {
                pl.data()[j] = 0;
            }
        }

        Ciphertext tmp1, tmp2;
        s1 = chrono::high_resolution_clock::now();
        evaluator.multiply_plain(ct1, pl, tmp1);
        evaluator.multiply_plain(ct2, pl, tmp2);
        e1 = chrono::high_resolution_clock::now();
        // cout << "        after multi plain two time: " << chrono::duration_cast<chrono::microseconds>(e1 - s1).count() << endl;
        evaluator.transform_from_ntt_inplace(tmp1);
        evaluator.transform_from_ntt_inplace(tmp2);

        s1 = chrono::high_resolution_clock::now();
        if (i == 0){
            evaluator.multiply(tmp1, tmp2, final_ct);
        } else {
            Ciphertext tmp;
            evaluator.multiply(tmp1, tmp2, tmp);
            evaluator.add_inplace(final_ct, tmp);
        } 
        e1 = chrono::high_resolution_clock::now();
        // cout  << "        after multi...time: " << chrono::duration_cast<chrono::microseconds>(e1 - s1).count() << endl;
    }

    evaluator.relinearize_inplace(final_ct, relin_keys);
    return final_ct;
}


inline
Ciphertext EvalAddMany_inpace_modImprove_extract_load(const int start_index, const SEALContext& context,
                                                      SecretKey& sk, int ct_size) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    chrono::high_resolution_clock::time_point s1, e1;
    uint64_t loading = 0;
    while (ct_size != 1) {
        for(int i = 0; i < ct_size/2; i++) {
            Ciphertext ct1, ct2;
            s1 = chrono::high_resolution_clock::now();
            loadCiphertext(context, ct1, start_index+i);
            loadCiphertext(context, ct2, start_index+ct_size/2+i);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            evaluator.add_inplace(ct1, ct2);
            s1 = chrono::high_resolution_clock::now();
            saveCiphertext(ct1, start_index+i);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
        }
        if (ct_size%2 == 0)
            ct_size = ct_size/2;
        else { // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            // ciphertexts[ct_size/2] = ciphertexts[ct_size-1];
            Ciphertext tmp;
            s1 = chrono::high_resolution_clock::now();
            loadCiphertext(context, tmp, start_index+ct_size-1);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            s1 = chrono::high_resolution_clock::now();
            saveCiphertext(tmp, start_index+ct_size/2);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            ct_size = ct_size/2+1;
        }
        counter++;
    }

    Ciphertext res;
    s1 = chrono::high_resolution_clock::now();
    loadCiphertext(context, res, start_index);
    e1 = chrono::high_resolution_clock::now();
    loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();

    loading_time += loading;
    // cout << "loading time: " << loading << endl;
    return res;
}


inline
Ciphertext EvalMultMany_inpace_modImprove_extract_load(const int start_index, const RelinKeys &relin_keys,
                                                       const SEALContext& context, SecretKey& sk, int ct_size) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    chrono::high_resolution_clock::time_point s1, e1;
    uint64_t loading = 0;
    while (ct_size != 1) {
        cout << "Level: " << ct_size << " , counter: " << counter << endl;
        for(int i = 0; i < ct_size/2; i++) {
            // cout << "   " << start_index+i << ", " << start_index+ct_size/2+i << endl;
            Ciphertext ct1, ct2;
            s1 = chrono::high_resolution_clock::now();
            loadCiphertext(context, ct1, start_index+i);
            loadCiphertext(context, ct2, start_index+ct_size/2+i);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            // s1 = chrono::high_resolution_clock::now();
            ct1 = extract_and_multiply(relin_keys, context, ct1, ct2);
            // e1 = chrono::high_resolution_clock::now();
            // cout << "---------------------------: " << chrono::duration_cast<chrono::microseconds>(e1 - s1).count() << endl;
            // evaluator.multiply_inplace(ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            // evaluator.relinearize_inplace(ciphertexts[i], relin_keys);
            if(counter & 1) {
                if (i == 0 ) cout << "  mod switched.\n";
                evaluator.mod_switch_to_next_inplace(ct1);
            }
            s1 = chrono::high_resolution_clock::now();
            saveCiphertext(ct1, start_index+i);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
        }
        if (ct_size%2 == 0)
            ct_size = ct_size/2;
        else { // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            // ciphertexts[ct_size/2] = ciphertexts[ct_size-1];
            Ciphertext tmp;
            s1 = chrono::high_resolution_clock::now();
            loadCiphertext(context, tmp, start_index+ct_size-1);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(tmp);
            }
            s1 = chrono::high_resolution_clock::now();
            saveCiphertext(tmp, start_index+ct_size/2);
            e1 = chrono::high_resolution_clock::now();
            loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();
            ct_size = ct_size/2+1;
        }
        counter++;
    }

    Ciphertext res;
    s1 = chrono::high_resolution_clock::now();
    loadCiphertext(context, res, start_index);
    e1 = chrono::high_resolution_clock::now();
    loading += chrono::duration_cast<chrono::microseconds>(e1 - s1).count();

    loading_time += loading;
    // cout << "loading time: " << loading << endl;
    return res;
}


inline
Ciphertext EvalMultMany_inpace_modImprove_extract_iterator(vector<Ciphertext>::iterator ciphertexts_it,
                                                           const RelinKeys &relin_keys, const SEALContext& context,
                                                           SecretKey& sk, int ct_size) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ct_size != 1){
        for(int i = 0; i < ct_size/2; i++){
            *(ciphertexts_it+i) = extract_and_multiply(relin_keys, context, *(ciphertexts_it+1), 
                                                       *(ciphertexts_it+ct_size/2+i));

            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(*(ciphertexts_it+i));
            }
        }
        if(ct_size%2 == 0)
            ct_size = ct_size/2;
        else{ // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            *(ciphertexts_it+ct_size/2) = *(ciphertexts_it+ct_size-1);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(*(ciphertexts_it+ct_size/2));
            }
            ct_size = ct_size/2+1;
        }
        counter += 1;
    }

    Ciphertext res = *ciphertexts_it;
    return res;
}


inline
Ciphertext EvalAddMany_inpace_modImprove_extract_multi_core(vector<Ciphertext> ciphertexts, const SEALContext& context,
                                                            SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
            evaluator.add_inplace(ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}


inline
Ciphertext EvalMultMany_inpace_modImprove_extract_multi_core(vector<Ciphertext> ciphertexts, const RelinKeys &relin_keys,
                                                             const SEALContext& context, SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
            ciphertexts[i] = extract_and_multiply_multi_core(relin_keys, context, ciphertexts[i],
                                                             ciphertexts[ciphertexts.size()/2+i]);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[i]);
            }
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[ciphertexts.size()/2]);
            }
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}


inline
Ciphertext EvalMultMany_inpace_modImprove_extract(vector<Ciphertext> ciphertexts, const RelinKeys &relin_keys,
                                                  const SEALContext& context, SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
            ciphertexts[i] = extract_and_multiply(relin_keys, context, ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[i]);
            }
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[ciphertexts.size()/2]);
            }
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}

inline
Ciphertext EvalMultMany_inpace_modImprove(vector<Ciphertext>& input, const RelinKeys &relin_keys,
                                          const SEALContext& context, SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    vector<Ciphertext> ciphertexts(input.size());
    for (int i = 0; i < (int) input.size(); i++) {
        ciphertexts[i] = input[i];
    }

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
            evaluator.multiply_inplace(ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            evaluator.relinearize_inplace(ciphertexts[i], relin_keys);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[i]);
            }
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[ciphertexts.size()/2]);
            }
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}


inline
long power(long x, long y, long m)
{
    if (y == 0)
        return 1;
    long p = power(x, y / 2, m) % m;
    p = (p * p) % m;
 
    return (y % 2 == 0) ? p : (x * p) % m;
}

inline
long modInverse(long a, long m)
{
    return power(a, m - 2, m);
}


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb,
                                                             const int n = 1024, const int p = prime_p,
                                                             const uint64_t big_prime = 1152921504589938689) {
    vector<regevCiphertext> results(ring_dim);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = std::make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());
    uniform_int_distribution<uint32_t> dist(0, 100);

    for (int cnt = 0; cnt < ring_dim; cnt++) {
        results[cnt].a = NativeVector(n);
        int ind = 0;
        for (int i = cnt; i >= 0 && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((long double) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = temp < 0 ? p + temp : temp;

            ind++;
        }

        for (int i = ring_dim-1; i > ring_dim - n + cnt && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((long double) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = -temp < 0 ? p-temp : -temp;

            ind++;
        }

        float temp_f = ((float) rlwe_ct.data(0)[cnt]) * ((float) p) / ((long double) big_prime);
        uint32_t decimal = temp_f - ((int) temp_f) * 100;
        float rounding = dist(engine) < decimal ? 1 : 0;

        long temp = ((int) (temp_f + rounding)) % p;
        results[cnt].b = temp % ((int) p);
    }

    return results;
}



inline void multiply_power_of_X(EncryptionParameters& enc_param, const Ciphertext &encrypted, Ciphertext &destination,
                                uint32_t index) {

    auto coeff_mod_count = enc_param.coeff_modulus().size() - 1;
    auto coeff_count = enc_param.poly_modulus_degree();
    auto encrypted_count = encrypted.size();

    destination = encrypted;

    for (int i = 0; i < (int) encrypted_count; i++) {
        for (int j = 0; j < (int) coeff_mod_count; j++) {
            negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count),
                                           coeff_count, index,
                                           enc_param.coeff_modulus()[j],
                                           destination.data(i) + (j * coeff_count));
        }
    }
}

inline vector<Ciphertext> subExpand(const SecretKey& sk, const SEALContext& context, EncryptionParameters& enc_param,
                                    const Ciphertext &encrypted, uint32_t m, const GaloisKeys& galkey,
                                    int first_expansion_size, int t = 65537) {

    Evaluator evaluator(context);
    Plaintext two("2");
    Decryptor decryptor(context, sk);

    int logFirst = ceil(log2(first_expansion_size));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = 0; i < logFirst; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        uint32_t index = (static_cast<uint64_t>(index_raw) * galois_elts[i]) % (m << 1);

        for (uint32_t a = 0; a < temp.size(); a++) {

            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
            // Let later children reuse this parent's buffer in the SEAL pool.
            temp[a].release();
        }

        temp.swap(newtemp);
    }

    temp.resize(first_expansion_size);
    return temp;
}

inline vector<Ciphertext> expand(const SEALContext& context, EncryptionParameters& enc_param, const Ciphertext &encrypted,
                                 uint32_t m, const GaloisKeys& galkey, int stepSize = 1, int t = 65537) {

    Evaluator evaluator(context);
    Plaintext two("2");

    int first_expansion_size = m / stepSize;
    int logFirst = ceil(log2(first_expansion_size));
    int logm = ceil(log2(m));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = logFirst; i < logm - 1; i++) {
        // cout << i << endl;
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (m << 1);

        for (int a = 0; a < (int) temp.size(); a++) {

            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);

            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);

            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }

        temp = newtemp;
    }

    // cout << "now?" << endl;

    // Last step of the loop
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (m << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (m << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) { // corner case.
            evaluator.multiply_plain(temp[a], two, newtemp[a]); // plain multiplication by 2.
        } else {
            evaluator.apply_galois(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + stepSize;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}



// for a tree with m leaf node, m >> stepSize, we first expand it to a subtree with m / stepSize leaf node
// (i.e., this subtree is the top of the whole tree)
// and then for each leaf node in this subtree, expand it into a small subtree with stepSize leaf node
inline vector<Ciphertext> expand_standalone(const SEALContext& context, EncryptionParameters& enc_param,
                                 const SecretKey& sk, const Ciphertext &encrypted, uint32_t m,
                                 const GaloisKeys& galkey, int t = 65537) {

    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    Plaintext two("2");
    int logm = ceil(log2(m));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = 0; i < logm - 1; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (m << 1);

        for (int a = 0; a < (int) temp.size(); a++) {
            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);

            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);

            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }

        temp = newtemp;
    }

    // Last step of the loop
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (m << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (m << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) { // corner case.
            evaluator.multiply_plain(temp[a], two, newtemp[a]); // plain multiplication by 2.
        } else {
            evaluator.apply_galois(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);
        
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + m;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}


vector<vector<uint64_t>> generateMatrixU_transpose(int n, const int q = 65537) {
    vector<vector<uint64_t>> U(n,  vector<uint64_t>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                U[i][j] = (uint64_t) power_seal(primitive_root, j, q);
            } else if (i == n/2) {
                U[i][j] = (uint64_t) modInverse_seal(U[0][j], q);
            } else {
                U[i][j] = (uint64_t) power_seal(U[i-1][j], 3, q);
            }
        }
    }
    return U;
}


Ciphertext coeffToSlot_WOPreprocess(const SEALContext& context, Ciphertext& input_ct, const GaloisKeys& gal_keys,
                                     const int degree, const uint64_t q = 65537, const uint64_t scalar = 1,
                                     const int numcores = 8) {

    int sq_rt = sqrt(degree/2);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    
    vector<Ciphertext> input_sqrt_list(2*sq_rt);

    chrono::high_resolution_clock::time_point time_start, time_end, s,e;
    uint64_t total_U = 0, total_multi = 0;

    time_start = chrono::high_resolution_clock::now();

    Ciphertext input_ct_copy(input_ct);
	evaluator.rotate_columns_inplace(input_ct_copy, gal_keys);

    input_sqrt_list[0] = input_ct;
	input_sqrt_list[sq_rt] = input_ct_copy;

    for (int c = 1; c < sq_rt; c++) {
        evaluator.rotate_rows(input_sqrt_list[c-1], sq_rt, gal_keys, input_sqrt_list[c]);
        evaluator.rotate_rows(input_sqrt_list[c-1+sq_rt], sq_rt, gal_keys, input_sqrt_list[c+sq_rt]);
    }

    NTL::SetNumThreads(numcores);
    NTL_EXEC_RANGE(sq_rt, first, last);
    for (int i = first; i < last; i++) {
        evaluator.transform_to_ntt_inplace(input_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(input_sqrt_list[i+sq_rt]);
    }
    NTL_EXEC_RANGE_END;

    time_end = chrono::high_resolution_clock::now();
    cout << "   ** [TIME] ** prepare rotated ciphertext for CoeffToSlot: "\
         << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    
    // time_start = chrono::high_resolution_clock::now();
    // vector<vector<uint64_t>> U = generateMatrixU_transpose(degree, q);
    // cout << "After generating U...\n";
    // vector<vector<uint64_t>> U_inverse = generateInverse_vander(U, q);
    // cout << "After generating the inverse of matrix U...\n";
    // saveDoubleVector(U_inverse);
    
    time_start = chrono::high_resolution_clock::now();
    vector<vector<uint64_t>> U_inverse(degree, vector<uint64_t>(degree));
    loadDoubleVector(U_inverse);

    time_end = chrono::high_resolution_clock::now();
    total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    vector<Ciphertext> result(sq_rt);
    int count = 0;

    NTL::SetNumThreads(numcores);
    NTL_EXEC_RANGE(sq_rt, first, last);
    for (int iter  = first; iter < last; iter++) {
            for (int j = 0; j < (int) input_sqrt_list.size(); j++) {

                s = chrono::high_resolution_clock::now();
                time_start = chrono::high_resolution_clock::now();
                vector<uint64_t> U_tmp(degree);
                for (int i = 0; i < degree; i++) {
                    int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                    row_index = i < degree/2 ? row_index : row_index + degree/2;
                    int col_index = (i + j*sq_rt) % (degree/2);
                    if (j < (int) input_sqrt_list.size() / 2) { // first half
                        col_index = i < degree/2 ? col_index : col_index + degree/2;
                    } else {
                        col_index = i < degree/2 ? col_index + degree/2 : col_index;
                    }
                    U_tmp[i] = ((uint64_t) (U_inverse[row_index][col_index] * scalar)) % q;
                }

                Plaintext U_plain;
                batch_encoder.encode(U_tmp, U_plain);
                evaluator.transform_to_ntt_inplace(U_plain, input_sqrt_list[j].parms_id());

                time_end = chrono::high_resolution_clock::now();
                U_time_multi_core += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

                count++;
                if (j == 0) {
                    evaluator.multiply_plain(input_sqrt_list[j], U_plain, result[iter]);
                } else {
                    Ciphertext temp;
                    evaluator.multiply_plain(input_sqrt_list[j], U_plain, temp);
                    evaluator.add_inplace(result[iter], temp);
                }
                e = chrono::high_resolution_clock::now();
                total_multi += chrono::duration_cast<chrono::microseconds>(e - s).count();
            }
    }
    NTL_EXEC_RANGE_END;

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    U_time += total_U;
    cout << "   ** [TIME] ** multiplication time: " << total_multi << endl;

    return result[0];
}

vector<Ciphertext> coeffToSlot_WOPreprocess_batch(SecretKey& sk, SEALContext& context, vector<Ciphertext>& input_ct, const GaloisKeys& gal_keys,
                                                  const int degree, const int batch_size, const uint64_t q = 65537, const uint64_t scalar = 1,
                                                  const int numcores = 8) {

    // assume each ct has batch_size info, we need iter ciphertexts to pack them together
    int iter = ceil((double) (input_ct.size() * batch_size) / (double) degree);
    int pack_in_half_ct = floor( (double) (degree/2) / (double) (batch_size));
    int pack_in_one_ct = 2*pack_in_half_ct;

    vector<Ciphertext> results(iter);

    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    BatchEncoder batch_encoder(context);

    int core_share = input_ct.size() / numcores;

    vector<vector<Ciphertext>> core_results(results.size(), vector<Ciphertext>(numcores));
    NTL::SetNumThreads(numcores);

    NTL_EXEC_RANGE(input_ct.size(), first, last);

    for (int i = first; i < last; i++) {
        Plaintext rotator;
        rotator.resize(degree);
        rotator.parms_id() = parms_id_zero;
        for (int j = 0; j < degree; j++) {
            rotator.data()[j] = 0;
        }
        int index_ct = floor((double) i / (double) pack_in_one_ct); // the index of result ciphertext
        int rot_row_offset = batch_size * (i % pack_in_half_ct); // slot on row
        bool rot_col_offset = (i % pack_in_one_ct) >= pack_in_half_ct ? true : false; // if should rotate the first half and second half
        int rot_offset = rot_row_offset + degree/2 * rot_col_offset; // encoded started from this slot
        rotator.data()[rot_offset] = 1;

        if (i % core_share == 0) {
            evaluator.multiply_plain(input_ct[i], rotator, core_results[index_ct][i/core_share]);
        } else {
            Ciphertext extracted_input;
            evaluator.multiply_plain(input_ct[i], rotator, extracted_input);
            evaluator.add_inplace(core_results[index_ct][i/core_share], extracted_input);
        }
    }

    NTL_EXEC_RANGE_END;

    for (int iter = 0; iter < (int) results.size(); iter++) {
        for (int i = 0; i < numcores; i++) {
            if (i == 0) {
                results[iter] = core_results[iter][i];
            } else {
                evaluator.add_inplace(results[iter], core_results[iter][i]);
            }
        }
    }


    for (int i = 0; i < (int) results.size(); i++) {
        results[i] = coeffToSlot_WOPreprocess(context, results[i], gal_keys, degree, q, scalar, numcores);
    }
    print_ct_to_vec(results[0], context, sk, degree);

    // vector<Ciphertext> unpack_results(2 * results.size()); // extract one-out-of-two consecutive packed selection vectors, for poly eval
    // Plaintext extractor_pl;
    // vector<uint64_t> extractor(degree);
    // for (int i = 0; i < (int) results.size(); i++) {
    //     for (int j = 0; j < degree/2; j++) {
    //         extractor[j] = (j / batch_size) % 2 == 0 ? 1 : 0; // the 0-th batch, 2-nd batch, 4-th batch, ....
    //         extractor[j+degree/2] = extractor[j];
    //     }
    //     batch_encoder.encode(extractor, extractor_pl);
    //     evaluator.multiply_plain(results[i], extractor_pl, unpack_results[2*i]);

    //     for (int j = 0; j < degree/2; j++) {
    //         extractor[j] = (j / batch_size) % 2 == 1 ? 1 : 0; // the 1-st batch, 3-rd batch, 5-th batch, ....
    //         extractor[j+degree/2] = extractor[j];
    //     }
    //     batch_encoder.encode(extractor, extractor_pl);
    //     evaluator.multiply_plain(results[i], extractor_pl, unpack_results[2*i + 1]);
    // }

    return results;
}



Ciphertext slotToCoeff_WOPrepreocess(const SEALContext& context, Ciphertext& input_ct, const GaloisKeys& gal_keys,
                                     const int degree, const uint64_t q = 65537, const uint64_t scalar = 1,
                                     const int numcores = 8, const int group_size = 32768, const int ring_dim = 32768) {

    int sq_rt = sqrt(degree/2);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_U = 0, total = 0;

    time_start = chrono::high_resolution_clock::now();
    vector<Ciphertext> input_sqrt_list(2*sq_rt);

    Ciphertext input_ct_copy(input_ct);
	evaluator.rotate_columns_inplace(input_ct_copy, gal_keys);

    input_sqrt_list[0] = input_ct;
	input_sqrt_list[sq_rt] = input_ct_copy;

    for (int c = 1; c < sq_rt; c++) {
        evaluator.rotate_rows(input_sqrt_list[c-1], sq_rt, gal_keys, input_sqrt_list[c]);
        evaluator.rotate_rows(input_sqrt_list[c-1+sq_rt], sq_rt, gal_keys, input_sqrt_list[c+sq_rt]);
    }
    for (int c = 0; c < sq_rt; c++) {
        evaluator.transform_to_ntt_inplace(input_sqrt_list[c]);
        evaluator.transform_to_ntt_inplace(input_sqrt_list[c+sq_rt]);
    }
    time_end = chrono::high_resolution_clock::now();

    total += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    
    time_start = chrono::high_resolution_clock::now();
    vector<vector<uint64_t>> U = generateMatrixU_transpose(degree, q);
    time_end = chrono::high_resolution_clock::now();
    total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    preprocess_time += total_U;

    struct PhaseTimes {
        uint64_t prepare_us;
        uint64_t evaluate_us;
    };
    vector<PhaseTimes> times_by_partition(numcores);
    const auto total_s = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> result(sq_rt);

    int core_share = sq_rt / numcores;
    NTL_EXEC_RANGE(numcores, first, last);
    for (int c  = first; c < last; c++) {
        PhaseTimes times{};
        for (int iter = c * core_share; iter < (c+1) * core_share; iter++) {
            for (int j = 0; j < (int) input_sqrt_list.size(); j++) {

                auto start = chrono::high_resolution_clock::now();
                vector<uint64_t> U_tmp(degree);
                for (int i = 0; i < degree; i++) {
                    int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                    row_index = i < degree/2 ? row_index : row_index + degree/2;
                    int col_index = (i + j*sq_rt) % (degree/2);
                    if (j < (int) input_sqrt_list.size() / 2) { // first half
                        col_index = i < degree/2 ? col_index : col_index + degree/2;
                    } else {
                        col_index = i < degree/2 ? col_index + degree/2 : col_index;
                    }
                    U_tmp[i] = ((uint64_t) (U[row_index][col_index] * scalar)) % q;
                }

                
                Plaintext U_plain;
                batch_encoder.encode(U_tmp, U_plain);
                evaluator.transform_to_ntt_inplace(U_plain, input_sqrt_list[j].parms_id());
                auto end = chrono::high_resolution_clock::now();
                times.prepare_us += chrono::duration_cast<chrono::microseconds>(end - start).count();

                start = chrono::high_resolution_clock::now();
                if (j == 0) {
                    evaluator.multiply_plain(input_sqrt_list[j], U_plain, result[iter]);
                } else {
                    Ciphertext temp;
                    evaluator.multiply_plain(input_sqrt_list[j], U_plain, temp);
                    evaluator.add_inplace(result[iter], temp);
                }
                end = chrono::high_resolution_clock::now();
                times.evaluate_us += chrono::duration_cast<chrono::microseconds>(end - start).count();
            }
        }
        times_by_partition[c] = times;
    }
    NTL_EXEC_RANGE_END;

    const auto total_e = chrono::high_resolution_clock::now();
    PhaseTimes sums{};
    for (const auto& times : times_by_partition) {
        sums.prepare_us += times.prepare_us;
        sums.evaluate_us += times.evaluate_us;
    }
    const auto elapsed_us = chrono::duration_cast<chrono::microseconds>(total_e - total_s).count();
    const auto measured_us = sums.prepare_us + sums.evaluate_us;
    total += elapsed_us * sums.evaluate_us / measured_us;
    preprocess_time += elapsed_us * sums.prepare_us / measured_us;

    time_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }
    time_end = chrono::high_resolution_clock::now();
    total += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // cout << "       second slotToCoeff time: " << total << endl;
    
    total_time += total / (ring_dim/group_size); // amortized

    return result[0];
}


Ciphertext calculateDegree(const SEALContext& context, const RelinKeys &relin_keys, Ciphertext& input, map<int, bool> modDownIndices, int degree) {
	Evaluator evaluator(context);

	vector<Ciphertext> output(degree);
	output[0] = input;
	vector<int> calculated(degree, 0), numMod(degree, 0);
	calculated[0] = 1;

	Ciphertext res, base;

	auto toCalculate = degree;
	int resdeg = 0;
	int basedeg = 1;
	base = input;
	while(toCalculate > 0){
		if(toCalculate & 1){
			toCalculate -= 1;
			resdeg += basedeg;
			if(calculated[resdeg-1] != 0){
				res = output[resdeg - 1];
			} else {
				if(resdeg == basedeg){
					res = base; // should've never be used as base should have made calculated[basedeg-1]
				} else {
					numMod[resdeg-1] = numMod[basedeg-1];

					evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
					evaluator.multiply_inplace(res, base);
					evaluator.relinearize_inplace(res, relin_keys);
					if(modDownIndices.count(resdeg) && !modDownIndices[resdeg]) {
						modDownIndices[resdeg] = true;
						evaluator.mod_switch_to_next_inplace(res);
						numMod[resdeg-1]+=1;
					}
				}
				output[resdeg-1] = res;
				calculated[resdeg-1] += 1;
			}
		} else {
			toCalculate /= 2;
			basedeg *= 2;
			if(calculated[basedeg-1] != 0){
				base = output[basedeg - 1];
			} else {
				numMod[basedeg-1] = numMod[basedeg/2-1];
				evaluator.square_inplace(base);
				evaluator.relinearize_inplace(base, relin_keys);
				while(modDownIndices.count(basedeg) && !modDownIndices[basedeg]) {
					modDownIndices[basedeg] = true;
					evaluator.mod_switch_to_next_inplace(base);
					numMod[basedeg-1]+=1;
				}
				output[basedeg-1] = base;
				calculated[basedeg-1] += 1;
			}
		}
	}

	return output[output.size()-1];
}


Ciphertext raisePowerToPrime(const SEALContext& context, const RelinKeys &relin_keys, Ciphertext& input, map<int, bool> modDownIndices_1,
							 map<int, bool> modDownIndices_2, int degree_1, int degree_2, int prime = 65537) {

	Ciphertext tmp = calculateDegree(context, relin_keys, input, modDownIndices_1, degree_1);
	tmp = calculateDegree(context, relin_keys, tmp, modDownIndices_2, degree_2);

	return tmp;
}


vector<vector<uint64_t>> generateEvaluationMatrix(int start_index, const int ring_dim, const int q = 65537) {
    vector<vector<uint64_t>> evaluation(ring_dim, vector<uint64_t>(ring_dim));

    for (int i = 0; i < ring_dim; i++) {
        for (int j = 0; j < ring_dim; j++) {
            if (i == 0) {
                evaluation[i][j] = 1;
            } else if (j == 0) {
                evaluation[i][j] = 1;
            } else {
                evaluation[i][j] = power_seal(start_index + i+1, j, q);
            }
        }
    }

    return evaluation;
}


Ciphertext evaluatePolynomial(const SEALContext context, const Ciphertext poly_ct, const GaloisKeys& gal_keys,
                              const int ring_dim, const int evaluate_index, const int numcores = 1,
                              const int q = 65537) {
    int sq_rt = sqrt(ring_dim/2);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    
    vector<Ciphertext> input_sqrt_list(2*sq_rt);

    Ciphertext poly_ct_copy(poly_ct);
    evaluator.rotate_columns_inplace(poly_ct_copy, gal_keys);

    input_sqrt_list[0] = poly_ct;
	input_sqrt_list[sq_rt] = poly_ct_copy;

    chrono::high_resolution_clock::time_point ss, ee;

    ss = chrono::high_resolution_clock::now();
    for (int c = 1; c < sq_rt; c++) {
        evaluator.rotate_rows(input_sqrt_list[c-1], sq_rt, gal_keys, input_sqrt_list[c]);
        evaluator.rotate_rows(input_sqrt_list[c-1+sq_rt], sq_rt, gal_keys, input_sqrt_list[c+sq_rt]);
    }
    for (int c = 0; c < sq_rt; c++) {
        evaluator.transform_to_ntt_inplace(input_sqrt_list[c]);
        evaluator.transform_to_ntt_inplace(input_sqrt_list[c+sq_rt]);
    }
    ee = chrono::high_resolution_clock::now();
    cout << "###### prepare time: " << chrono::duration_cast<chrono::microseconds>(ee - ss).count() << endl;

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_U = 0, total_multi = 0;

    vector<Ciphertext> result(sq_rt);

    int count = 0;

    for (int eval_iter = 0; eval_iter < evaluate_index/ring_dim + 1; eval_iter++) {

        time_start = chrono::high_resolution_clock::now();
        vector<vector<uint64_t>> evaluation_mat = generateEvaluationMatrix(ring_dim * eval_iter, ring_dim);
        time_end = chrono::high_resolution_clock::now();
        total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

        int core_share = sq_rt / numcores;
        NTL_EXEC_RANGE(numcores, first, last);
        for (int c  = first; c < last; c++) {
            for (int iter = c * core_share; iter < (c+1) * core_share; iter++) {
                for (int j = 0; j < (int) input_sqrt_list.size(); j++) {
                    ss = chrono::high_resolution_clock::now();
                    time_start = chrono::high_resolution_clock::now();
                    vector<uint64_t> eval_tmp(ring_dim);
                    for (int i = 0; i < ring_dim; i++) {
                        int row_index = (i-iter) % (ring_dim/2) < 0 ? (i-iter) % (ring_dim/2) + ring_dim/2 : (i-iter) % (ring_dim/2);
                        row_index = i < ring_dim/2 ? row_index : row_index + ring_dim/2;
                        int col_index = (i + j*sq_rt) % (ring_dim/2);
                        if (j < (int) input_sqrt_list.size() / 2) { // first half
                            col_index = i < ring_dim/2 ? col_index : col_index + ring_dim/2;
                        } else {
                            col_index = i < ring_dim/2 ? col_index + ring_dim/2 : col_index;
                        }
                        eval_tmp[i] = ((uint64_t) (evaluation_mat[0][0])) % q;
                    }

                    Plaintext U_plain;
                    batch_encoder.encode(eval_tmp, U_plain);
                    evaluator.transform_to_ntt_inplace(U_plain, input_sqrt_list[j].parms_id());

                    time_end = chrono::high_resolution_clock::now();
                    U_time_multi_core += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

		            count++;
                    if (eval_iter == 0 && j == 0) {
                        evaluator.multiply_plain(input_sqrt_list[j], U_plain, result[iter]);
                    } else {
                        Ciphertext temp;
                        evaluator.multiply_plain(input_sqrt_list[j], U_plain, temp);
                        evaluator.add_inplace(result[iter], temp);
                    }

                    ee = chrono::high_resolution_clock::now();
                    total_multi += chrono::duration_cast<chrono::microseconds>(ee - ss).count();
                }
            }
        }
        NTL_EXEC_RANGE_END;
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    U_time += total_U;
    cout << "###### " << count << " multiplications.\n";
    cout << "   ** [TIME] ** multiplication time: " << total_multi << ", " << total_multi / count << endl;


    return result[0];
}


vector<vector<Ciphertext>> evaluatePolynomial_batch(SecretKey& sk, SEALContext context, const vector<Ciphertext> poly_ct, const GaloisKeys& gal_keys,
                                            const int ring_dim, const int evaluatePoly_batch_size, const int evaluatePoly_party_size,
                                            const int evaluatePoly_degree, const int numcores = 8, const int q = 65537) {
                                                
    int sq_batch = sqrt(evaluatePoly_batch_size); // 20 for batch_size = sbar = 400
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    int pack_in_half_ct = floor( (double) (ring_dim/2) / (double) (2*evaluatePoly_batch_size));

    vector<uint64_t> extractor(ring_dim);
    vector<Plaintext> extractor_ntt(sq_batch);
    for (int i = 0; i < sq_batch; i++) {
        for (int j = 0; j < ring_dim/2; j++) {
            if (j >= 2*evaluatePoly_batch_size * pack_in_half_ct) {
                extractor[j] = 0;
                extractor[j + ring_dim/2] = 0;
            } else if ((j / evaluatePoly_batch_size) % 2 == 1) {
                extractor[j] = 0;
                extractor[j + ring_dim/2] = 0;
            } else {
                extractor[j] = ((j % evaluatePoly_batch_size) / sq_batch) == i ? 1 : 0;
                extractor[j + ring_dim/2] =  extractor[j];
            }
        }
        batch_encoder.encode(extractor, extractor_ntt[i]);
        evaluator.transform_to_ntt_inplace(extractor_ntt[i], poly_ct[0].parms_id());
    }

    vector<vector<Ciphertext>> results(poly_ct.size(), vector<Ciphertext>(164));
    
    int flip = 0; // indicate whether the non-zero batch is in the 2nd or 1st allocate batch_size slots
    for (int poly_ind = 0; poly_ind < (int) poly_ct.size(); poly_ind++) {

        // for each ct, same as normal evaluation w/o batching, we optimize the ntt and rotation and multiplication, prepare sq of ntt first
        vector<Ciphertext> poly_ntt(sq_batch);
        Ciphertext tmp;
        evaluator.multiply_plain(poly_ct[poly_ind], extractor_ntt[0], tmp);
        evaluator.rotate_rows(tmp, -evaluatePoly_batch_size, gal_keys, poly_ntt[0]);
        evaluator.add_inplace(poly_ntt[0], poly_ct[poly_ind]);
        for (int i = 1; i < sq_batch; i++) {
            evaluator.multiply_plain(poly_ct[poly_ind], extractor_ntt[sq_batch-i], tmp);
            evaluator.rotate_rows_inplace(tmp, (sq_batch - i)*sq_batch, gal_keys);
            evaluator.rotate_rows(poly_ntt[i-1], -sq_batch, gal_keys, poly_ntt[i]);
            evaluator.add_inplace(poly_ntt[i], tmp);
        }
        
        for (int i = 0; i < sq_batch; i++) {
            evaluator.transform_to_ntt_inplace(poly_ntt[i]);
        }

        chrono::high_resolution_clock::time_point s1, e1;
        s1 = chrono::high_resolution_clock::now();
        // 65537 values, 164 batches, each 400 values, the last batch take 400-63 values
        NTL_EXEC_RANGE(164, first, last);
        for (int c = first; c < last; c++) {
            vector<uint64_t> indices(evaluatePoly_batch_size); // prepare the indices values
            for (int i = 0; i < evaluatePoly_batch_size; i++) {
                indices[i] = c*evaluatePoly_batch_size + i;
                indices[i] = indices[i] >= 65537 ? 0 : indices[i];
            }

            vector<Ciphertext> result_one_batch(sq_batch);
            for (int i = 0; i < sq_batch; i++) {
                vector<uint64_t> tmp_vec(ring_dim);

                for (int j = 0; j < sq_batch; j++) {
                    Plaintext indice_pl;
                    for (int k = 0; k < ring_dim/2; k++) {
                        // fall in the allocated rotation slots, always zero
                        if ((k < i) || (((k-i) / evaluatePoly_batch_size) % 2 != flip) || (k >= 2*evaluatePoly_batch_size * pack_in_half_ct)) {
                            tmp_vec[k] = 0;
                        } else {
                            int ind = (k-i) % evaluatePoly_batch_size;
                            int power = ind - j*sq_batch < 0 ? ind - j*sq_batch + evaluatePoly_batch_size : ind - j*sq_batch;
                            power = (power + i) % evaluatePoly_batch_size;
                            tmp_vec[k] = power_seal(indices[ind], power, q);
                        }
                        tmp_vec[k+ring_dim/2] = tmp_vec[k];
                    }
                    batch_encoder.encode(tmp_vec, indice_pl);
                    evaluator.transform_to_ntt_inplace(indice_pl, poly_ntt[0].parms_id());

                    if (j == 0) {
                        evaluator.multiply_plain(poly_ntt[j], indice_pl, result_one_batch[i]);
                    } else {
                        Ciphertext tmp_ct;
                        evaluator.multiply_plain(poly_ntt[j], indice_pl, tmp_ct);
                        evaluator.add_inplace(result_one_batch[i], tmp_ct);
                    }
                }
            }

            for (int i = 0; i < sq_batch; i++) {
                evaluator.transform_from_ntt_inplace(result_one_batch[i]);
            }

            for (int i = sq_batch-1; i >=0; i--) {
                if (i == sq_batch-1) {
                    results[poly_ind][c] = result_one_batch[i];
                } else {
                    evaluator.rotate_rows_inplace(results[poly_ind][c], 1, gal_keys); 
                    evaluator.add_inplace(results[poly_ind][c], result_one_batch[i]);
                }
            }
        }
        NTL_EXEC_RANGE_END;
        e1 = chrono::high_resolution_clock::now();
        cout << "Total eval time: " << chrono::duration_cast<chrono::microseconds>(e1 - s1).count() << endl;

        flip = 1 - flip;
    }

    return results;
}


void expandFirstToAll(SecretKey& sk, SEALContext context, Ciphertext& input, const GaloisKeys& gal_keys,
                      const int ring_dim) {

    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    // vector<Plaintext> extractor(log2(ring_dim));
    // vector<uint64_t> vvv(ring_dim);
    // NTL_EXEC_RANGE(log2(ring_dim), first, last);
    // for (int i = first; i < last; i++) {
    //     for (int j = 0; j < ring_dim; j++) {
    //         vvv[j] = 0;
    //     }
    //     int gap = ring_dim / (1<<i);
    //     for (int j = 0; j < (1<<i); j++) {
    //         vvv[gap * j] = 1;
    //     }
    //     batch_encoder.encode(vvv, extractor[i]);
    // }
    // NTL_EXEC_RANGE_END;

    // int counter = 0;
    
    // Plaintext extractor;
    for (int i = 0; i < log2(ring_dim); i++) {
        // for (int j = 0; j < ring_dim; j++) {
        //     vvv[j] = 0;
        // }
        int gap = ring_dim / (1<<i);
        // for (int j = 0; j < (1<<i); j++) {
        //     vvv[gap * j] = 1;
        // }
        // batch_encoder.encode(vvv, extractor);
        Ciphertext extracted(input);
        // evaluator.transform_to_ntt_inplace(extractor[i], input.parms_id());
        // evaluator.multiply_plain(input, extractor[i], extracted);
        // // evaluator.transform_from_ntt_inplace(extracted);

        int rotate_ind = gap / 2;
        if (rotate_ind == ring_dim / 2) {
            evaluator.rotate_columns_inplace(extracted, gal_keys);
        } else {
            evaluator.rotate_rows_inplace(extracted, -rotate_ind, gal_keys);
        }
        evaluator.add_inplace(input, extracted); 

        // counter++;
        // if (counter == 3) {
        //     counter = 0;
        //     evaluator.mod_switch_to_next_inplace(input);
        // }
    }
}


// Ciphertext sumUpEvaluation_ToPartyRange(SEALContext& context, vector<Ciphertext>& evaluated, const int group_size,
//                                         GaloisKeys& gal_keys, const int ring_dim = 32768, const int eval_size = 65536) {
//     Ciphertext result;

//     Evaluator evaluator(context);

//     int iter = log2(ring_dim / group_size); // need to fold iter times to map eval -> group

//     evaluator.add(evaluated[0], evaluated[1], result);

//     for (int i = 0; i < iter; i++) {
//         int rot_ind = ring_dim / (1<<(i+1)); // i=0, fold half, rot=16384; i=1, rot=8192
//         Ciphertext tmp;
//         if (rot_ind == ring_dim/2) {
//             evaluator.rotate_columns(result, gal_keys, tmp);
//         } else {
//             evaluator.rotate_rows(result, rot_ind, gal_keys, tmp);
//         }
//         evaluator.add_inplace(result, tmp);
//     }

//     return result;
// }


void sumUpEvaluation_ToPartyRange(SEALContext& context, Ciphertext& evaluated, const int group_size,
                                        GaloisKeys& gal_keys, const int ring_dim = 32768, const int eval_size = 65536) {
    // Ciphertext result;

    Evaluator evaluator(context);

    int iter = log2(ring_dim / group_size); // need to fold iter times to map eval -> group

    // evaluator.add(evaluated[0], evaluated[1], result);

    for (int i = 0; i < iter; i++) {
        int rot_ind = ring_dim / (1<<(i+1)); // i=0, fold half, rot=16384; i=1, rot=8192
        Ciphertext tmp;
        if (rot_ind == ring_dim/2) {
            evaluator.rotate_columns(evaluated, gal_keys, tmp);
        } else {
            evaluator.rotate_rows(evaluated, rot_ind, gal_keys, tmp);
        }
        evaluator.add_inplace(evaluated, tmp);
    }

    // return result;
}

