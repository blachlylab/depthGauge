module csv;

import std.stdio;
import std.algorithm:splitter,joiner;
import std.conv:to;
import std.array:array,appender,join;
import std.algorithm.iteration:each,map;
import std.algorithm.sorting:sort;
import std.range;
import std.traits:ReturnType;
import bio.std.hts.sam.header:SamHeader;
import bio.std.hts.bam.reader:BamReader;
import std.path;
struct Table {
    string[] header;
    string[] samples;
    Record[] records;
    uint[string] contigs;
    SamHeader samheader;

    File f;
    string delim;
    ReturnType!createMatrix matrix;
    this(string filename,string prefix,int startSamples){
        this(filename,startSamples,prefix,"\t");
    }
    this(string filename,int startSamples,string prefix, string delim){
        f=File(filename);
        this.delim=delim;
        parse(startSamples,prefix);
    }
    void parse(int startSamples,string prefix){
        auto lines=f.byLineCopy();
        //set header
        header=lines.front.splitter(delim).array;
        //create samples
        samples=header[startSamples..$];
        auto bam=new BamReader(buildPath(prefix,samples[0]~".bam"));
        samheader=bam.header();
        //debug writeln(header);
        //debug writeln(samples);
        lines.popFront;
        //create records
        foreach(line;lines){
            auto split=line.splitter(delim);
            auto rec=split.take(startSamples).array;
            records~=Record(rec);
        }
        matrix=createMatrix();
    }
    void write(File f){
        f.writeln(join(header,delim));
        foreach(i,rec;enumerate(records.sort)){
            f.writeln(join([samheader.getSequence(rec.chr).name,rec.pos.to!(string)]~rec.extra~matrix[i][].map!(x=>x.to!(string)).array,delim));
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
    this(string[] line,ref SamHeader samheader){
        chr=samheader.getSequenceIndex(line.front);
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
