// ==========================================================================
//               LaRAgu - Lagrangian Relaxation Aligner GU
// ==========================================================================
// Copyright (c) 2015-2016, Gianvito Urgese
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Gianvito Urgese nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL GIANVITO URGESE OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
//         Orazio Scicolone
// ==========================================================================

#ifndef STOPWATCH_H
#define STOPWATCH_H


#include <chrono>

template<typename TimeT = std::chrono::microseconds,
        typename ClockT=std::chrono::high_resolution_clock,
        typename DurationT=double>
class Stopwatch
{
private:
    std::chrono::time_point<ClockT> _start, _end;
    std::string _label;
public:
    Stopwatch(std::string str) {
        start();
        this->_label= str;
    }
    void start() {
        _start = _end = ClockT::now();
    }
    DurationT stop() {
        _end = ClockT::now(); return elapsed();
    }
    DurationT elapsed() {
        auto delta = std::chrono::duration_cast<TimeT>(_end-_start);
        return delta.count();
    }
//    void print(std::ofstream & out){
//        this->stop();
//        out << this->_label << " (Sec) = " << std::to_string(this->elapsed()/1000000) << std::endl;
//    }
    void print(std::ostream & out){
        this->stop();
        out << this->_label << " (Sec) = " << std::to_string(this->elapsed()/1000000) << std::endl;
    }
};


#endif //STOPWATCH_H
