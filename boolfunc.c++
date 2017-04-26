#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "weight_table.h"
#include <time.h>
#include "weight_table_v2.h"


#define BASE unsigned long long

using namespace std;

class BF 
{
private:
		BASE * pointer;
		int n_args;
public:
		BF();
    BF(const BF &);
    BF(const char *);
    BF(int type, int args);// type: 0 - const 0; 1 - const 1; 2 - random;
    ~BF();
    void set(int x, int pos);// set x in location pos;
    int get(int pos);
    int args();
    BF& operator = (const BF &);
    bool operator == (const BF &);
    void print(const char * name = "F");
    void print_all(const char * name = "F");
    int weight();
    int weight_v2();
    int BASE_weight(int );
    int weight_v3();
    int weight_v4();
    void print_anf(const char * name = "F");
    BF mebius();
    int deg();
    void PUA(int *);
    int cor_imm();
    int nonlinearity();
    BF NAP();
    void auto_corr(long long *);
    int SAC_1();
    int SAC();
    int PC();
};

///////////////////////////////////////////////////////////////////////////////////////////////

    void BF::set(int type, int pos)
{
    if((pos < 0) || (pos >= (1<<n_args))) cout<<"\nBF::set(int, int): wrong position.\n";
    else
    {
          // cout<<"\npos = "<<pos<<"; pos/64 = "<<pos/64<<"; (1<<(pos%64)) = "<<((BASE)1<<(pos%64)); 
        //cout<<"\nbefore = "<<(BASE)pointer[pos/64]<<"  after = "<<(BASE)(pointer[pos/64] | ((BASE)1<<(pos%64)));
        if(type == '1')  pointer[pos/64] = pointer[pos/64] | ((BASE)1<<(pos%64));
        if(type == 1)    pointer[pos/64] = pointer[pos/64] | ((BASE)1<<(pos%64));
    }
}

/*-------------------------------------------------------------------------------------------*/

    int BF::get(int pos)
{
    if((pos < 0) || (pos >= (1<<n_args))) cout<<"\nBF::get( int): wrong position.\n";
    else
    {
    return (pointer[pos/64] & ((BASE)1<<(pos%64))) >> (pos%64);
    }
}

/*-------------------------------------------------------------------------------------------*/

    int BF::args()
{
    return n_args;
}

/*-------------------------------------------------------------------------------------------*/

    BF::BF()
{
    pointer=NULL;
    n_args=0;
}    

/*-------------------------------------------------------------------------------------------*/

    BF::BF(const BF & source)
{
    n_args=source.n_args;
    if(source.n_args<=6) 
    { 
        pointer = new BASE[1];
        pointer[0]=source.pointer[0]; 
    }
    else
    {
        pointer = new BASE[1<<(source.n_args-6)];
        for(int i=0; i<(1<<(source.n_args-6)); i++) {pointer[i]=source.pointer[i];}
    }
}

/*-------------------------------------------------------------------------------------------*/

    BF::~BF()
{
    if (pointer != NULL) delete [] pointer;
}

/*-------------------------------------------------------------------------------------------*/

    BF::BF(int type, int args)
{
    if (type == 0) 
    {
        n_args=args;
        if (args<=6) {pointer = new BASE[1]; pointer[0]=(BASE)0;}
        else 
        {
            pointer = new BASE[1<<(args-6)];
            for(int i=0; i<(1<<(args-6)); i++) pointer[i]=(BASE)0;
        } 
    }
    if (type == 1) 
    {
        n_args=args;
        if (args<=6) 
        {
            pointer = new BASE[1]; 
            if(args<6) pointer[0]=((BASE)1<<(1<<args))-1;
            else pointer[0]=0xFFFFFFFFFFFFFFFF;
        }
        else 
        {
            pointer = new BASE[1<<(args-6)]; 
            for(int i=0; i<(1<<(args-6)); i++) pointer[i]=0xFFFFFFFFFFFFFFFF;
        }
    }
    if (type == 2)
    {
        n_args=args;
        if (args<=6) 
        {
            pointer = new BASE[1]; 
            if(args<6) pointer[0]=(((BASE)(rand()-rand())<<32) + (BASE)(rand()-rand()))&(((BASE)1<<(1<<args))-1);
            else pointer[0]=((BASE)(rand()-rand())<<32) + (BASE)(rand()-rand());
        }
        else 
        {
            pointer = new BASE[1<<(args-6)]; 
            for(int i=0; i<(1<<(args-6)); i++) 
                pointer[i]=((BASE)(rand()-rand())<<32)+(BASE)(rand()-rand());
        }
    }
}

/*-------------------------------------------------------------------------------------------*/
    
    BF::BF(const char * source)
{
    if (source == NULL) {n_args=0; pointer=NULL;}
    else 
    {
        int size = strlen(source);
        n_args = 0;
        if ((size & (size-1)) == 0)
        {
            while((size&1)!=1) {n_args++; size=size>>1;}
            if (n_args<=6) pointer = new BASE[1];
            else pointer = new BASE[1<<(n_args-6)];
            for(int i=0; i<(1<<n_args); i++) set(source[i],i);
        }
        else {cout<<"\nWrong lenght of function!\n"; pointer=NULL; }
    }
}

/*-------------------------------------------------------------------------------------------*/

    BF& BF::operator = (const BF &source)
{
    n_args=source.n_args;
    if(n_args<=6) {pointer = new BASE[1]; pointer[0] = source.pointer[0];}
    else 
    {
        pointer = new BASE[1<<(n_args-6)];
        for(int i=0; i<(1<<(n_args-6)); i++) pointer[i] = source.pointer[i];
    }
    return *this;
}

/*-------------------------------------------------------------------------------------------*/

    bool BF::operator == (const BF &source)
{
    if(n_args != source.n_args) return 0;
    if(n_args<=6) if(pointer[0] != source.pointer[0]) return 0;
    else return 1;
    for(int i=0; i<(1<<(n_args-6)); i++) if(pointer[i] != source.pointer[i]) return 0;
    return 1;
}

/*-------------------------------------------------------------------------------------------*/

    void BF::print(const char * name)
{
    //head
    int len=strlen(name)|1;
    cout<<endl<<"+";
    for(int i=0; i<n_args+5+len; i++) cout<<"-";
    cout<<"+"<<endl<<"| ";
    for(int i=1; i<=n_args; i++) cout<<i;    
    cout<<" | "<<name;
    if((strlen(name)&0x1)==0) cout<<" ";
    cout<<" |"<<endl<<"+";
    for(int i=0; i<n_args+2; i++) cout<<"-";
    cout<<"+";
    for(int i=0; i<len+2; i++) cout<<"-";
    cout<<"+"<<endl;
    //body

    for(int i=0; i<(1<<n_args); i++)
    {
        cout<<"| ";
        for(int j=n_args-1; j>=0; j--) cout<<((i&(1<<j))>>j);
        cout<<" | ";
        for(int j=0; j<(len>>1); j++) cout<<" ";
        cout<<get(i);
        for(int j=0; j<(len>>1); j++) cout<<" ";
        cout<<" |"<<endl;
    }
    //endl
    
    cout<<"+";
    for(int i=0; i<n_args+5+len; i++) cout<<"-";
    cout<<"+"<<endl;
}

/*-------------------------------------------------------------------------------------------*/

    void BF::print_all(const char * name)
{
    int * pua = new int[1<<n_args];
    long long * AC = new long long[1<<n_args];
    this->PUA(pua); this->auto_corr(AC);

    //head

    int len=strlen(name)|1;
    cout<<endl<<"+";
    for(int i=0; i<n_args+18+len; i++) cout<<"-";
    cout<<"+"<<endl<<"| ";
    for(int i=1; i<=n_args; i++) cout<<i;    
    cout<<" | "<<name;
    if((strlen(name)&0x1)==0) cout<<" ";
    cout<<" | PUA |  AC  |"<<endl<<"+";
    for(int i=0; i<n_args+2; i++) cout<<"-";
    cout<<"+";
    for(int i=0; i<len+2; i++) cout<<"-";
    cout<<"+-----+------+"<<endl;
    //body

    for(int i=0; i<(1<<n_args); i++)
    {
        cout<<"| ";
        for(int j=n_args-1; j>=0; j--) cout<<((i&(1<<j))>>j);
        cout<<" | ";
        for(int j=0; j<(len>>1); j++) cout<<" ";
        cout<<get(i);
        for(int j=0; j<(len>>1); j++) cout<<" ";
        cout<<" |";
        printf(" %3d | %4lld |\n", pua[i], AC[i]);
    }
    //endl
    
    cout<<"+";
    for(int i=0; i<n_args+18+len; i++) cout<<"-";
    cout<<"+"<<endl;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::weight_v3()
{
    if(n_args<=6) 
    {
        return W[pointer[0]&0xFF] + W[(pointer[0]>>8)&0xFF] + W[(pointer[0]>>16)&0xFF]+ W[(pointer[0]>>24)&0xFF] + W[(pointer[0]>>32)&0xFF] + W[(pointer[0]>>40)&0xFF] + W[(pointer[0]>>48)&0xFF] + W[(pointer[0]>>56)&0xFF]; 
    }
    else 
    {
        int tmp=0;
        for(int i=0; i<(1<<(n_args-6)); i++) 
        {
            tmp += W[pointer[i]&0xFF] + W[(pointer[i]>>8)&0xFF] + W[(pointer[i]>>16)&0xFF]+ W[(pointer[i]>>24)&0xFF] + W[(pointer[i]>>32)&0xFF] + W[(pointer[i]>>40)&0xFF] + W[(pointer[i]>>48)&0xFF] + W[(pointer[i]>>56)&0xFF];
        }
        return tmp;
    }
}

/*-------------------------------------------------------------------------------------------*/

    int BF::BASE_weight(int base_n)
{
    //    cout<<endl<<"BASE_weight["<<base_n<<"]"<<endl;
    BASE tmp = pointer[base_n];
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = tmp - ((tmp>>1) & 0x5555555555555555);
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = (tmp & 0x3333333333333333) + ((tmp>>2) & 0x3333333333333333);
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = (tmp + (tmp>>4)) & 0x0f0f0f0f0f0f0f0f;
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = tmp + (tmp>>8);
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = tmp + (tmp>>16);
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = tmp + (tmp>>32);
    //    for(int i=63; i>=0; i--) cout<<((tmp&((BASE)1<<i))>>i);
    //    cout<<endl;
    tmp = tmp & 0x7f;
    return (int)tmp;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::weight_v2()
{
    if(n_args<=6) return this->BASE_weight(0);
    else
    {
        int tmp=0;
        for(int i=0; i<(1<<(n_args-6)); i++) 
        {
            //cout<<"tmp_before = "<<tmp<<endl;
            //cout<<"BASE_w["<<i<<"]  = "<<this->BASE_weight(i)<<endl; 
            tmp = tmp + this->BASE_weight(i);
            //cout<<"tmp_after  = "<<tmp<<endl<<endl;
        }
        return tmp;
    }
}

/*-------------------------------------------------------------------------------------------*/

    int BF::weight()
{
    if(n_args<=8) return this->weight_v3();
    else
    {
        int V=0, i=0;
        int BASE_N = 1<<(n_args-6);
        while(i<BASE_N)    
        {
            int gr, j=i;
            BASE s=0;
            if((i+31) < BASE_N) gr = i+31;
            else gr = BASE_N;
            while(j<gr)
            {
                BASE x = pointer[j];
                x = x - ((x>>1) & 0x5555555555555555);
                x = (x & 0x3333333333333333) + ((x>>2) & 0x3333333333333333);
                x = (x + (x>>4)) & 0x0f0f0f0f0f0f0f0f;
                s = s + x; j++;
            }
            s = (s & 0x00ff00ff00ff00ff) + ((s>>8) & 0x00ff00ff00ff00ff);
            s = (s & 0x0000ffff0000ffff) + ((s>>16) & 0x0000ffff0000ffff);
            s = (s & 0x00000000ffffffff) + (s>>32);
            V = V + (int)s; i+=31;
        }
        return V;
    }
}

/*-------------------------------------------------------------------------------------------*/

    int BF::weight_v4()
{
    if(n_args<=6) 
    {
        return WW[pointer[0]&0xFFFF] + WW[(pointer[0]>>16)&0xFFFF] + WW[(pointer[0]>>32)&0xFFFF]+ WW[(pointer[0]>>48)&0xFFFF]; 
    }
    else 
    {
        int tmp=0;
        for(int i=0; i<(1<<(n_args-6)); i++) 
        {
            tmp +=  WW[pointer[i]&0xFFFF] + WW[(pointer[i]>>16)&0xFFFF] + WW[(pointer[i]>>32)&0xFFFF]+ WW[(pointer[i]>>48)&0xFFFF]; 
        }
        return tmp;
    }
}

/*-------------------------------------------------------------------------------------------*/

    BF BF::mebius()
{
    BF tmp(*this);
    if(n_args <=6) 
    {
        tmp.pointer[0] = tmp.pointer[0] ^ ((tmp.pointer[0]<< 1) & 0xaaaaaaaaaaaaaaaa);
        tmp.pointer[0] = tmp.pointer[0] ^ ((tmp.pointer[0]<< 2) & 0xcccccccccccccccc);
        tmp.pointer[0] = tmp.pointer[0] ^ ((tmp.pointer[0]<< 4) & 0xf0f0f0f0f0f0f0f0);
        tmp.pointer[0] = tmp.pointer[0] ^ ((tmp.pointer[0]<< 8) & 0xff00ff00ff00ff00);
        tmp.pointer[0] = tmp.pointer[0] ^ ((tmp.pointer[0]<<16) & 0xffff0000ffff0000);
        tmp.pointer[0] = tmp.pointer[0] ^ (tmp.pointer[0]<<32);
    }
    else
    {
        int k = 1<<(n_args-6), j = 2;
        for(int i=0; i<(1<<(n_args-6)); i++) 
        {
            tmp.pointer[i] = tmp.pointer[i] ^ ((tmp.pointer[i]<< 1) & 0xaaaaaaaaaaaaaaaa);
            tmp.pointer[i] = tmp.pointer[i] ^ ((tmp.pointer[i]<< 2) & 0xcccccccccccccccc);
            tmp.pointer[i] = tmp.pointer[i] ^ ((tmp.pointer[i]<< 4) & 0xf0f0f0f0f0f0f0f0);
            tmp.pointer[i] = tmp.pointer[i] ^ ((tmp.pointer[i]<< 8) & 0xff00ff00ff00ff00);
            tmp.pointer[i] = tmp.pointer[i] ^ ((tmp.pointer[i]<<16) & 0xffff0000ffff0000);
            tmp.pointer[i] = tmp.pointer[i] ^ (tmp.pointer[i]<<32);
        }
        while(k>1)
        {
            for(int q=0; q<(1<<(n_args-6)); q+=j)
                for(int i=0; i<(j>>1); i++) tmp.pointer[q+i+(j>>1)] ^= tmp.pointer[q+i];
            k = k>>1; j = j<<1;
        }
    }
    return tmp;
}

/*-------------------------------------------------------------------------------------------*/

    void BF::print_anf(const char * name)
{
    if(n_args<=10) {
    BF tmp = this->mebius();
    cout<<endl<<name<<" = ";
    int q=0;
    if(tmp.get(0)) {cout<<"1"; q++;}
    for(int i=1; i<(1<<tmp.n_args); i++)
    {
        if(tmp.get(i))
        {
            if(q > 0) cout<<" + ";
            for(int j=tmp.n_args-1, k=1; j>=0; j--, k++)
            {
                if((i&(1<<j)) > 0) cout<<"X"<<k;
            }
            q++;
        }
    }
    if(q == 0) cout<<"0";
    cout<<endl<<endl; }
    else cout<<endl<<"Too many args."<<endl;
}

/*-------------------------------------------------------------------------------------------*/

int BF::deg()
{
    BF tmp = this->mebius();
    if(tmp.get((1<<tmp.n_args)-1)) return tmp.n_args;
    /*int D=0, k=0;
    for(int i=(1<<tmp.n_args)-2; i>=0; i--) 
    {
    
        if(tmp.get(i))
        {
            int j=i;
            k=0;
            while(j>0) {j=j&(j-1); k++;}
            if(k > D) {
                    D = k;
            }
        }
    }
    return D;*/
    for(int k=n_args-1; k>=0; k--)
    {
        int b, d, c, w;
        d = (1<<k)-1;
        b = d<<(n_args-k);
        if(tmp.get(b)) return k;
        while(b!=d)
        {
            c = (b+1)&b;
            w = 0;
            int q=(c-1)^b;
            while(q>0) {q=q&(q-1); w++;}
            w = w-2;
            b = (((((b+1)^b)<<1)+1)<<w)^c;
            if(tmp.get(b)) return k;
        }
    }
}

/*-------------------------------------------------------------------------------------------*/

    void BF::PUA(int * pua)
{
    int * tmp = new int[1<<n_args];
    for(int i=0; i<(1<<n_args); i++) if(get(i)) pua[i] = -1; else pua[i] = 1;
    for(int i=0, j=1; i<n_args; i++, j=j<<1) 
    {
        for(int k=0; k<(1<<n_args); k=k+(j<<1)) 
        {
            for(int l=0; l<j;      l++) tmp[k+l]   = pua[k+l] + pua[k+l+j];
            for(int l=0; l<j;      l++) tmp[k+l+j] = pua[k+l] - pua[k+l+j];
            for(int l=0; l<(j<<1); l++) pua[k+l]   = tmp[k+l];
        }
    }
    delete []tmp;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::cor_imm()
{
    BF const_1(1,n_args), const_0(0,n_args);
    if(*this == const_1 || *this == const_0) return n_args;

    int * pua = new int[1<<n_args]; 
    this->PUA(pua);
    
    for(int k=1; k<=n_args-1; k++)
    {
        int b, d, c, w;
        d = (1<<k)-1;
        b = d<<(n_args-k);
        if(pua[b] != 0) return k-1;
        while(b!=d)
        {
          	c = (b+1)&b;
        	  w = 0;
        	  int q=(c-1)^b;
         	  while(q>0) {q=q&(q-1); w++;}
         	  w = w-2;
         	  b = (((((b+1)^b)<<1)+1)<<w)^c;
        	  if(pua[b] != 0) return k-1;
        }
    }
    delete []pua;
    return n_args-1;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::nonlinearity()
{
    int * pua = new int [1<<n_args];
    this->PUA(pua);
    int max=0;
    for(int i=0; i<(1<<n_args); i++)
    {
        if(abs(pua[i]) > max) max = abs(pua[i]);
    }
    delete []pua;
    return (1<<(n_args-1)) - (max>>1);
}    

/*-------------------------------------------------------------------------------------------*/

    BF BF::NAP()
{
    int * pua = new int [1<<n_args];
    this->PUA(pua);
    int max=0, max_ind;
    for(int i=0; i<(1<<n_args); i++)
    {
        if(abs(pua[i]) > max) {max = abs(pua[i]); max_ind = i;}
    }
    BF nap(0,n_args);
    for(int i=n_args; i>=0; i--) if((max_ind & (1<<i))) nap.set(1,(max_ind & (1<<i)));
    if(pua[max_ind] < 0) nap.set(1,0);
    delete []pua;
    return nap.mebius();
}

/*-------------------------------------------------------------------------------------------*/

    void BF::auto_corr(long long * AC)
{
    long long * tmp = new long long[1<<n_args];
    for(int i=0; i<(1<<n_args); i++) if(get(i)) AC[i] = -1; else AC[i] = 1;
    for(int i=0, j=1; i<n_args; i++, j=j<<1) 
    {
        for(int k=0; k<(1<<n_args); k=k+(j<<1)) 
        {
            for(int l=0; l<j;      l++) tmp[k+l]   = AC[k+l] + AC[k+l+j];
            for(int l=0; l<j;      l++) tmp[k+l+j] = AC[k+l] - AC[k+l+j];
            for(int l=0; l<(j<<1); l++) AC[k+l]    = tmp[k+l];
        }
    }
  
    for(int i=0; i<(1<<n_args); i++) AC[i] = AC[i] * AC[i];
    
    for(int i=0, j=1; i<n_args; i++, j=j<<1) 
    {
        for(int k=0; k<(1<<n_args); k=k+(j<<1)) 
        {
            for(int l=0; l<j;      l++) tmp[k+l]   = AC[k+l] + AC[k+l+j];
            for(int l=0; l<j;      l++) tmp[k+l+j] = AC[k+l] - AC[k+l+j];
            for(int l=0; l<(j<<1); l++) AC[k+l]    = tmp[k+l]>>1;
        }
    }
    delete []tmp;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::SAC_1()
{
    long long * AC = new long long[1<<n_args];
    this->auto_corr(AC);
    for(int i=0; i<n_args; i++) if(AC[1<<i] != 0) {delete []AC; return 0;}
    delete []AC;
    return 1;
}

/*-------------------------------------------------------------------------------------------*/
/*  trans_ind() usage: in_size = sizeof(base); if wrong combination of pos and val return -1 */

    int trans_ind(int in_size, int out_size, int pos, int val, int base)
{
    int ind[in_size - out_size], v[in_size - out_size];
    for(int i=0, j=0; i<in_size; i++) if(((pos & (1<<i))>>i) == 1) { ind[j] = i; j++;}
    for(int i=0; i<(in_size - out_size); i++) v[i] = (val & (1<<i))>>i; 
    for(int i=0; i<(in_size - out_size); i++)
        if(((base & (1<<ind[i]))>>ind[i]) != v[i]) return -1;
    int new_base = base, tmp_h, tmp_l;
    for(int i=0; i<(in_size - out_size); i++)
    {
        tmp_l = new_base & ((1<<(ind[i]-i)) - 1);
        tmp_h = (new_base & (0x7FFFFFFF ^ ((1<<((ind[i]-i) + 1)) - 1)))>>1;
        new_base = tmp_h + tmp_l;
    }
    return new_base;
}

/*-------------------------------------------------------------------------------------------*/

    int BF::SAC()
{
    if(this->SAC_1())
    {
         //   cout<<"\nSAC_1 check!";
        for(int i=1; i<=(n_args-2); i++)              //generating all combinations of i-args
        {
            int b, d, c, w, x;
            d = (1<<i)-1;
            b = d<<(n_args-i);
            // check all subfunctions with fixed i-args in b positions
            x=b;
            int ind[i];
            for(int k=0, l=0; k<n_args; k++) if(((b & (1<<k))>>k) == 1) { ind[l] = k; l++;}
            while(x>0)
            {
                BF tmp(0,n_args-i);                   //create subfunction
                for(int j=0; j<(1<<n_args); j++)
                {
                    if((j & b) == x) 
                    {
                        int new_base = j, tmp_h, tmp_l;
                        for(int k=0; k<i; k++)
                        {
                            tmp_l = new_base & ((1<<(ind[k]-k)) - 1);
                            tmp_h = (new_base & (0x7FFFFFFF ^ ((1<<((ind[k]-k) + 1)) - 1)))>>1;
                            new_base = tmp_h + tmp_l;
                        }
//        cout<<"\n###################### BEFORE ####################################\nj     = ";
  //      for(int y=16; y>=0; y--) cout<<((j&((BASE)1<<y))>>y);
   //     cout<<"      F(j)    = "<<get(j)<<"\nnew_j = ";
    //    for(int y=16; y>=0; y--) cout<<((new_base&((BASE)1<<y))>>y);
     //   cout<<"  subF(new_j) = "<<tmp.get(new_base)<<endl;
                        tmp.set(get(j), new_base);
     //   cout<<"\n\n           AFTER\nj     = ";
      //  for(int y=16; y>=0; y--) cout<<((j&((BASE)1<<y))>>y);
       // cout<<"      F(j)    = "<<get(j)<<"\nnew_j = ";
       // for(int y=16; y>=0; y--) cout<<((new_base&((BASE)1<<y))>>y);
       // cout<<"  subF(new_j) = "<<tmp.get(new_base)<<"\n##################################################################"<<endl;
                    }
                }
                //tmp.print_all();
                if(tmp.SAC_1() != 1) return i - 1;    //check subfunction  
                x = (x - 1) & b;
            }
            while(b!=d)
            {
                c = (b+1)&b;
                w = 0;
                int q=(c-1)^b;
                while(q>0) {q=q&(q-1); w++;}
                w = w-2;
                b = (((((b+1)^b)<<1)+1)<<w)^c;
                // check all subfunctions with fixed i-args in b positions
                x=b;
                int ind[i];
                for(int k=0, l=0; k<n_args; k++) if(((b & (1<<k))>>k) == 1) { ind[l] = k; l++;}
                while(x>0)
                {
                    BF tmp(0,n_args-i);                   //create subfunction
                    for(int j=0; j<(1<<n_args); j++)
                    {
                        if((j & b) == x) 
                        {
                            int new_base = j, tmp_h, tmp_l;
                            for(int k=0; k<i; k++)
                            {
                                tmp_l = new_base & ((1<<(ind[k]-k)) - 1);
                                tmp_h = (new_base & (0x7FFFFFFF ^ ((1<<((ind[k]-k) + 1)) - 1)))>>1;
                                new_base = tmp_h + tmp_l;
                            }
                            tmp.set(get(j), new_base);
                        }
                    }
                    if(tmp.SAC_1() != 1) return i - 1;    //check subfunction  
                    x = (x - 1) & b;
                }
            }
        }
        return n_args - 2;
    }
    else return 0;
}

/*-------------------------------------------------------------------------------------------*/
    int BF::PC()
{
    BF const_1(1,n_args), const_0(0,n_args);
    if(*this == const_1 || *this == const_0) return 0;

    long long * AC = new long long[1<<n_args]; 
    this->auto_corr(AC);
    
    for(int k=1; k<=n_args; k++)
    {
        int b, d, c, w;
        d = (1<<k)-1;
        b = d<<(n_args-k);
        if(AC[b] != 0) return k-1;
        while(b!=d)
        {
          	c = (b+1)&b;
        	  w = 0;
        	  int q=(c-1)^b;
         	  while(q>0) {q=q&(q-1); w++;}
         	  w = w-2;
         	  b = (((((b+1)^b)<<1)+1)<<w)^c;
        	  if(AC[b] != 0) return k-1;
        }
    }
    delete []AC;
    return n_args;
}

///////////////////////////////////////// TESTS ///////////////////////////////////////////////

    void test_1()
{
    for(int i=2; i<32; i++) 
    {
        BF tmp(2,i);
        cout<<"args = "<<i<<"; weight/2^args = "<<(float)tmp.weight()/(float)((BASE)1<<i)<<endl;
    }
}

/*-------------------------------------------------------------------------------------------*/

    void test_2()
{
    float t=0;
    for(int i=0; i<100; i++)
    {
        int args=2+rand()%30;
        BF tmp(2, args);
        t+=(float)tmp.weight()/(float)((BASE)1<<args);
    }
    cout<<endl<<"avg(w/2^n) = "<<t/(float)100<<endl;
}

/*-------------------------------------------------------------------------------------------*/

    void test_3(int args)
{
    clock_t t_tmp, t_v1=0, t_v2=0, t_v3=0, t_v4=0;
    for(int i=0; i<10; i++)
    {
        BF tmp(2, args);

        t_tmp=clock();
        tmp.weight();
        t_tmp=clock()-t_tmp;
        t_v1+=t_tmp;

        t_tmp=clock();
        tmp.weight_v2();
        t_tmp=clock()-t_tmp;
        t_v2+=t_tmp;

        t_tmp=clock();
        tmp.weight_v3();
        t_tmp=clock()-t_tmp;
        t_v3+=t_tmp;

        t_tmp=clock();
        tmp.weight_v4();
        t_tmp=clock()-t_tmp;
        t_v4+=t_tmp;
    }
    cout<<"\nAVG time of BF::weight()    on "<<args<<" args = "<<((float)t_v1/(float)10)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<"AVG time of BF::weight_v2() on "<<args<<" args = "<<((float)t_v2/(float)10)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<"AVG time of BF::weight_v3() on "<<args<<" args = "<<((float)t_v3/(float)10)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<"AVG time of BF::weight_v4() on "<<args<<" args = "<<((float)t_v4/(float)10)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<"time_2/time_1 = "<<(float)t_v2/(float)t_v1<<endl<<endl;
    cout<<"time_3/time_1 = "<<(float)t_v3/(float)t_v1<<endl<<endl;
    cout<<"time_4/time_1 = "<<(float)t_v4/(float)t_v1<<endl<<endl;
}

/*-------------------------------------------------------------------------------------------*/

    void test_4()
{
    for(int i=0; i<32; i++)
    {   
        int args=i;    
        BF tmp(2,args);
        int v, v2, v3, v4;
        v=tmp.weight();
        v2=tmp.weight_v2();
        v3=tmp.weight_v3();
        v4=tmp.weight_v4();
        if(v!=v2) cout<<endl<<"ERROR in v2 on args = "<<args<<endl;
        if(v!=v3) cout<<endl<<"ERROR in v3 on args = "<<args<<endl;
        if(v!=v4) cout<<endl<<"ERROR in v4 on args = "<<args<<endl;
    }
}

/*-------------------------------------------------------------------------------------------*/

    void test_SAC(int args)
{
    BF tmp(0, args);
    int b, d, c, w, k=2;
    d = (1<<k)-1;
    b = d<<(args-k);
    tmp.set(1,b);
    while(b!=d)
    {
      	c = (b+1)&b;
     	  w = 0;
     	  int q=(c-1)^b;
     	  while(q>0) {q=q&(q-1); w++;}
     	  w = w-2;
     	  b = (((((b+1)^b)<<1)+1)<<w)^c;
        tmp.set(1,b);
    }
  //  tmp.print();
    tmp = tmp.mebius();
    //tmp.print_all();
    int q;
    clock_t tmp_t;
    tmp_t = clock();
    q=tmp.SAC();
    tmp_t = clock() - tmp_t;
    cout<<"\nTime of BF::SAC() on "<<args<<" args = "<<((float)tmp_t)/CLOCKS_PER_SEC<<"sec\nSAC(F) = "<<q<<endl;
//    tmp.print_anf();
}

//////////////////////////////////////// MAIN /////////////////////////////////////////////////

int main()
{
    srand(time(NULL));
    BF f(2,31);
  /*  BF f("00010110011010000110100010000000");
   // BF f1("00010111");
    //f.print();
    f=f.mebius();
   // f=f1;
    f.print_all();
    f.print_anf();
    cout<<"deg(F) = "<<f.deg()<<endl<<endl;
    cout<<"Corr_imm(F) = "<<f.cor_imm()<<endl<<endl;
    cout<<"Nonlinearity(F) = "<<f.nonlinearity()<<endl;
    (f.NAP()).print_anf("NAP(F)");
    cout<<"PC(F) = "<<f.PC()<<endl<<endl;
    cout<<"SAC(F) = "<<f.SAC()<<endl<<endl;
*/
//    for(int i=6; i<32; i++) test_SAC(i);


  //  test_SAC(8);

    //    PUA print
    int *pua = new int[1<<f.args()];
    f.PUA(pua); cout<<"PUA(F) = ";
    for(int i=0; i<(1<<f.args()); i++) cout<<pua[i]<<" ";
    cout<<endl<<endl;
    delete []pua;

    //    AC print
    /*long long *AC = new long long[1<<f.args()];
    f.auto_corr(AC); cout<<"AC(F) = ";
    for(int i=0; i<(1<<f.args()); i++) cout<<AC[i]<<" ";
    cout<<endl<<endl;
    delete []AC;*/

    return 0;
}
