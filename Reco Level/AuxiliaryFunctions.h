#ifndef AUXILIARYFUNCTIONS_H_
#define AUXILIARYFUNCTIONS_H_

namespace aux
{
template <class T>
bool greaterPt(T a, T b) {
    return a.pt() > b.pt(); }


template <class T>
bool lesserPt(T a, T b) {
    return a.pt() < b.pt(); }

template <class T>
bool isSameParticle(T a, T b) {
    return a.index == b.index; }
}

#endif
