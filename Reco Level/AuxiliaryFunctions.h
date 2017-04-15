#ifndef AUXILIARYFUNCTIONS_H_
#define AUXILIARYFUNCTIONS_H_

#include "DataFormats/Math/interface/deltaR.h"

namespace aux
{
template <class T>
inline bool greaterPt(T a, T b) {
    return a.pt() > b.pt(); }


template <class T>
inline bool lesserPt(T a, T b) {
    return a.pt() < b.pt(); }

template <class T>
inline bool isSameParticle(T a, T b) {
    return a.index == b.index; }

template <class T1, class T2>
inline double calcDeltaR(T1 a, T2 b) {
    return deltaR(a.eta(), a.phi(), b.eta(), b.phi());}

template <class T>
inline bool lesserFirstElement(T a, T b) {
    return a.first < b.first; }
}
#endif
