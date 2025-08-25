function z() {
  g++ -std=c++17 -O2 -DLOCAL -D_GLIBCXX_DEBUG "$1.cpp" 
}
function r() {
  ./a.out < "in$1"
}