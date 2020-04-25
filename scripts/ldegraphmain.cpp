#include <iostream>
#include <vector>
#include <string>

#include "../src/ldegraphalg.h"
#include "../src/ldealg.h"

using namespace std;

int main() {
  int x;
  vector<int> l;
  vector<int> r;
  while (cin >> x, x) {
    l.push_back(x);
  }
  while (cin >> x, x) {
    r.push_back(x);
  }
  LDEGraphAlg alg(l, r);
  cout << alg.solve().size() << '\n';
  return 0;
}
