#include<bits/stdc++.h>
using namespace std;

typedef long long ll;

#define rep(i, n, args...) for(int i = 0, ##args; i < n; ++i)
#define repr(i, s, e, args...) for(int i = s, ##args; i < e; ++i)
#define repl(i, s, e, args...) for(ll i = s, ##args; i < e; i++)
#define erg(i, n, args...) for(int i = vtx[n], ##args; ~i; i = eg[i].nxt)
#define fulset(x, v) memset(x, v, sizeof(x))
#define fulcpy(v, u) memset(v, u, sizeof(u));
#define SZ(x,t) (x) * sizeof(t)
#define OP(o, args...) operator o (##args) const
#define MF_graph(_N, _M, args...)\
	const static int N = _M, M = _M;\
	struct star {int v, nxt; ##args} eg[M];\
	int vtx[N], n, ec;\
	void init(int _n) { n = _n, ec = 0, memset(vtx, 0, SZ(n+1, int)); }
#define MF_addeg(v, tv, args...) eg[ec] = {v, tv, ##args}, tv = ec

const int INF = 0x3f3f3f3f, MOD = 1e9+7;
