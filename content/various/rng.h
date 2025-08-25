
/**
 * Author: Unknown
 * Date: Unknown
 * License: Unknown
 * Source: Unknown
 * Description: Random number generator using Mersenne Twister 64-bit algorithm
 * Time: O(1) initialization
 * Status: tested
 */
mt19937_64 rng((u32) chrono::steady_clock::now().time_since_epoch().count());
