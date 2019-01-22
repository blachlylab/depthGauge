module csv;

import std.stdio;
import std.algorithm:splitter,joiner;
import std.algorithm.sorting:sort;
import std.conv:to;
import std.array:array,appender,join;
import std.algorithm.iteration:each,map;
import std.range;
import std.traits:ReturnType;
import std.string:fromStringz;
import dhtslib.sam:SAMFile;

import std.path;

struct Table {
    string[] header;
    string[] samples;
    Record[] records;
    uint[string] contigs;
    SAMFile * sam;

    File f;
    string delim;
    ReturnType!createMatrix matrix;
    this(string filename,int startSamples){
        this(filename,startSamples,"\t");
    }
    this(string filename,int startSamples, string delim){
        f=File(filename);
        this.delim=delim;
        parseSamples(startSamples);
    }
    void parseSamples(int startSamples){
        auto lines=f.byLineCopy();
        //set header
        header=lines.front.splitter(delim).array;
        //create samples
        samples=header[startSamples..$];
    }
    void parseRecords(SAMFile * sam,int startSamples){
        auto lines=f.byLineCopy();
        this.sam=sam;
        //debug writeln(header);
        //debug writeln(samples);
        lines.popFront;
        //create records
        foreach(line;lines){
            auto split=line.splitter(delim);
            auto rec=split.take(startSamples).array;
            records~=Record(rec,sam);
        }
        matrix=createMatrix();
    }

    void write(File f){
        f.writeln(join(header,delim));
        foreach(i,rec;enumerate(records.sort)){
            f.writeln(join([fromStringz(sam.header.target_name[rec.chr]).idup,rec.pos.to!(string)]~rec.extra~matrix[i][].map!(x=>x.to!(string)).array,delim));
        }
    }
    auto createMatrix(){
        auto buf = new ulong[records.length * samples.length];
        return buf.chunks(samples.length);
    }
}

struct Sample{
    string name;
    this(string name){
        this.name=name;
    }
}

struct Record {
    int chr;
    //0-based
    uint pos;
    string[] extra;
    this(string[] line,SAMFile * sam){
        chr=sam.target_id(line.front);
        line.popFront;
        //convert from 1-based to 0-based
        pos=line.front.to!uint-1;
        line.popFront;
        extra=line.array;
    }
    int opCmp(const ref Record other) const nothrow{
        if(this.chr > other.chr){return 1;}
        if(this.chr <other.chr){return -1;}
        if(this.pos > other.pos){return 1;}
        if(this.pos < other.pos){return -1;}
        return 0;
    }
}

unittest{
    string[] test=["chr7".dup,"23483".dup,"A".dup,"T".dup,"V".dup,"odd".dup];
    auto r=Record(test);
    assert(r.chr=="chr7");
    assert(r.pos==23483);
    assert(r.extra==["A","T","V","odd"]);
}

//void main(string[] args){
//
//}
