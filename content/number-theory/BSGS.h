/**
 * Author: Bjorn Martinsson
 * Date: 2020-06-03
 * License: CC0
 * Source: own work
 * Description: Returns the smallest $N$ s.t. $X^N = Y \pmod P$, or
 * $-1$ if no such $N$ exists. modLog(X,1,P) can be used to
 * calculate the order of $X$.
 * Time: $O(\sqrt m)$
 * Status: tested for all 0 <= a,x < 500 and 0 < m < 500.
 *
 * Details: This algorithm uses the baby-step giant-step method to
 * find (i,j) such that a^(n * i) = b * a^j (mod m), where n > sqrt(m)
 * and 0 < i, j <= n. If a and m are coprime then a^j has a modular
 * inverse, which means that a^(i * n - j) = b (mod m$).
 *
 * However this particular implementation of baby-step giant-step works even
 * without assuming a and m are coprime, using the following idea:
 *
 * Assume p^x is a prime divisor of m. Then we have 3 cases
 *	 1. b is divisible by p^x
 *	 2. b is divisible only by some p^y, 0<y<x
 *	 3. b is not divisible by p
 * The important thing to note is that in case 2, modLog(a,b,m) (if
 * it exists) cannot be > sqrt(m), (technically it cannot be >= log2(m)).
 * So once all exponenents of a that are <= sqrt(m) has been checked, you
 * cannot have case 2. Case 2 is the only tricky case.
 *
 * So the modification allowing for non-coprime input involves checking all
 * exponents of a that are <= n, and then handling the non-tricky cases by
 * a simple gcd(a^n,m) == gcd(b,m) check.
 */
#pragma once

i64 BSGS(i64 X, i64 Y, i64 P) {
  X %= P, Y %= P;
  assert(gcd(X, P) == 1);
  const i64 B = sqrtl(P) + 1;
  unordered_map<i64, int> mp;
  for (i64 i = 0, cur = Y; i <= B; i++) {
    mp[cur] = i;
    cur = cur * X % P;
  }
  i64 step = 1;
  rep (i, 1, B) step = step * X % P;
  for (i64 p = 1, cur = 1; p <= B; p++) {
    cur = cur * step % P;
    if (mp.contains(cur)) {
      return B * p - mp[cur];
    }
  }
  return -1;
}

i64 exBSGS(i64 X, i64 Y, i64 P) {
  X %= P, Y %= P;
  i64 D = 1;
  int add = 0;
  while (true) {
    const auto g = gcd(X, P);
    if (g == 1) break;
    if (D == Y) return add;
    if (Y % g != 0) return -1;
    Y /= g;
    P /= g;
    add++;
    D = D * (X / g) % P;
  }
  const i64 B = sqrtl(P) + 1;
  unordered_map<i64, int> mp;
  for (i64 i = 0, cur = Y; i <= B; i++) {
    mp[cur] = i;
    cur = cur * X % P;
  }
  i64 step = 1;
  rep (i, 1, B) step = step * X % P;
  for (i64 p = 1, cur = D; p <= B; p++) {
    cur = cur * step % P;
    if (mp.contains(cur)) {
      return B * p - mp[cur] + add;
    }
  }
  return -1;
}