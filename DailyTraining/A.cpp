#include <bits/stdc++.h>

//#define int long long
#define endl '\n'
using namespace std;
using ll = long long;

template<typename A, typename B>
inline ostream &operator<<(ostream &out, const pair <A, B> &p) {
    return out << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
inline ostream &operator<<(ostream &out, const vector <T> &a) {
    out << "[";
    for (int i = 0; i < a.size(); i++) {
        if (i) out << ',';
        out << ' ' << a[i];
    }
    return out << " ]";
}

template<typename T>
inline ostream &operator<<(ostream &out, const set <T> &a) { return out << vector<T>(all(a)); }

const int N = 2e5 + 10, mod = 998244353;


signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    freopen("D:\\mypile\\acm\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acm\\out.txt", "w", stdout);

    int n; cin >> n;
    vector<int> a(n);
    vector<int> cnt(200010);
    int maxv = 0;
    for (int i = 0; i < n; i ++ ) {
        cin >> a[i];
        cnt[a[i]] ++ ;
        maxv = max(maxv, a[i]);
    }
    vector<vector<int>> dp(maxv + 1, vector<int>(4));

    for (int i = 1; i <= maxv; i ++ ) {
        for (int j = 0; j < 4; j ++ ) {
            dp[i][j] += dp[i - 1][j];
            if (cnt[i]) {
                dp[i][j] += dp[i - 1][j - 1] * cnt[i];
            }
        }
    }
    return 0;
}