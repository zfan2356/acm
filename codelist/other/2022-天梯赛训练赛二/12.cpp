//
// Created by onglu on 2022/3/9.
//

#include <bits/stdc++.h>

#define all(a) a.begin(),a.end()
#define rall(a) a.rbegin(),a.rend()

#define endl '\n'
#define lson (rt << 1)
#define rson (rt << 1 | 1)
#define Mid ((l + r) / 2)
//#define int long long
using namespace std;
const int N = 2e6 + 1009;
//const int N = 2e5 + 1009;
//const int N = 5009;
//const int N = 309
int n, m, a[N];
vector<int> ver[N];
int fa[N], d[N];
int fid(int x) {
    if(d[x] != 0) return d[x];
    if(fa[x] == 0) return 1;
    else return d[x] = fid(fa[x]) + 1;

}
void work() {
    cin >> n;
    for(int i = 1; i <= n; i++) {
        int k;
        cin >> k;
        for(int j = 0; j < k; j++) {
            int y;
            cin >> y;
            fa[y] = i;
        }
    }
    int p, ans = -1;
    for(int i = 1; i <= n; i++) {
        if(fid(i) > ans) {
            p = i;
            ans = fid(i);
        }
    }
    cout << p << endl;
}

signed main() {
#ifdef LOCAL
    freopen("C:\\Users\\onglu\\CLionProjects\\acm\\data.in", "r", stdin);
    freopen("C:\\Users\\onglu\\CLionProjects\\acm\\data.out", "w", stdout);
#endif
    ios::sync_with_stdio(false);
    cin.tie(0);
    work();
    return 0;
}