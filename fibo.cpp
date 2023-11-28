#include<iostream>
using namespace std;

void fibo(int order){
   
   int a_0 = 0, a_1 = 1, a_rec = 0;
   
   cout<<" "<<a_1;
   
   for(int i = 0; i < order-1; i++){
   
      a_rec = a_0 + a_1;
      a_0 = a_1;
      a_1 = a_rec;
      
      cout<<" "<<a_rec;
   }
   
   cout<<endl;
}


int main(int argc, char **argv){

   int order;

   cin>>order;
   fibo(order);

   return 0;
}
