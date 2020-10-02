//
// Created by Robert Stolz on 6/27/17.
//

#include "gene.h"
#include <thread>

//constructors
using namespace std;

void Gene::parse_header(){
    unsigned long pos = 0;
    string name, remaining;
    //extract gene name
    pos = header.find('=');
    name = header.substr(0,pos);
     //need to further parse the name
    pos = header.find(' ');
    name = name.substr(1,pos);
    remaining = header.substr(pos+1, header.length());
    pos = name.find(' ');
    name = name.substr(0,pos);
    gene_name = name;
    //extract chromosome name
    pos = remaining.find('=');
    remaining = remaining.substr(pos+1, remaining.length());
    pos = remaining.find(':');
    position.chromosome = remaining.substr(0,pos);
    //for (auto & c: position.chromosome) c = toupper(c); //C++11 string toUpper
    std::transform(position.chromosome.begin(), position.chromosome.end(), position.chromosome.begin(), ::tolower);
    std::transform(position.chromosome.end()-1, position.chromosome.end(), position.chromosome.end()-1, ::toupper);
    remaining = remaining.substr(pos+1, remaining.length());
    //extract the start and stop locations
    pos = remaining.find("-");
    position.start_pos = stol(remaining.substr(0,pos));
    remaining = remaining.substr(pos+1, remaining.length());
    pos = remaining.find(" ");
    position.end_pos = stol(remaining.substr(0,pos));
    remaining = remaining.substr(pos+1, remaining.length());
    pos = remaining.find("STRAND=");
    remaining = remaining.substr(pos+7, remaining.length());
    pos = remaining.find(" ");
    position.strand = remaining.substr(0,pos);

    return;
}
//constructors and destructors
Gene::Gene(){

}

Gene::~Gene(){
    clear_structures();
}

//getters and setters
string Gene::getName(){
    return gene_name;
}

const string &Gene::getHeader() const {
    return header;
}

void Gene::setHeader(const string &header) {
    Gene::header = header;
}

const Loci &Gene::getPosition() const {
    return position;
}

void Gene::setPosition(const Loci &position) {
    Gene::position = position;
}

const vector<char, allocator<char>> &Gene::getSequence() const {
    return sequence;
}

vector<Structure>& Gene::getStructures(){
    return rloop_structures;
}

float Gene::compute_GC_skew(){

    if (sequence.size()){
        throw EmptyGeneException();
    }
    float Gs = 0.f;
    float Cs = 0.f;
    for (int i=0; i < sequence.size(); i++){
        if (sequence[i] == 'G'){
            Gs += 1.f;
        }
        else if (sequence[i] == 'C'){
            Cs += 1.f;
        }
    }
    return (Gs-Cs)/(Gs+Cs);
}

bool Gene::read_gene(ifstream& fastafile) { //need to test
    //initialize variables
    char c,p;
    //check that fstream is open and not eof, throw exceptions otherwise
    if (!fastafile.is_open()) {
        throw UnexpectedClosedFileException("Gene::read_gene");
    } else if (fastafile.eof()) {
        throw UnexpectedEOFException();
    }
    while (fastafile.get(c)) {
        //read the next character
        c = toupper(c);
        //if the character is the start of a header line
        if (c == '>') {
            //read until the end of the header line
            while (c != '\n' && c != '\r') {
                header.push_back(c);
                c = toupper(fastafile.get());
            }
            parse_header();
        }
        else if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
            //save to sequence vector
            sequence.push_back(c);
        }
            //else if encountered some white space
        else if (c == '\n' || c == ' ' || c == '\t' || c== '\r'){
            p = fastafile.peek();
            if (p == '>'){
                windower.set_sequence(sequence);
                return false;
            }
            /*else if (p == EOF) {
                return true;
            }*/
            continue;
        }
            //else the charicter is unrecognized
        else {
            throw InvalidSequenceDataException(c);
        }
    }
    windower.set_sequence(sequence);
    return true;
    //unexpected EOF
    //throw InvalidSequenceDataException();
}

void Gene::print_gene(){
    std::cout << header << std::endl;
    for (std::vector<char>::iterator it = sequence.begin(); it<sequence.end(); ++it){
        cout << *it;
    }
    std::cout << std::endl;
}


void Gene::threaded_compute_structures(int block_start, int block_end, int number_threads){
// block start and end define what section this thread is responsible for
// calculating
    // for i in range(block start, end) iterate over block of vector
        *i // pointer to the structure
    //calculate the first structure to get us going
    // then enter the loop to dynamically calculate the rest in theory would
    // not need to check because the blocks should define the dynamic
    // programming boundry

    //std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
    // may need to make a sequence iterator that is offset by the block start
    // and stop so that could mean calculating the "reach" of the block
    // looking at the end position of the first and last structures in the block



    // need to create threads first
    int number_of_blocks = rloop_structures.back().block_id;
    // number of blocks that need to be distributed amoung the threads
    // ideally in a continous fashion

    // want to create vectors of pointers or just have them access regions of
    // that are defined by the blocks

    // need to figure out what the start and stop is for each
    // interactor from the rloop structure vector that gets passed to the
    // compute in each thread I think would work best
    // then everything would stay in place no worries about ordering because
    // the threads are just modifying the values of the objects stored in
    // the rloop structure vector


    std::vector<thread> thread_pool;
    for (int i; i < number_threads; i++){

        // need to figure out the chucks here and what chuncks will
        // be added to each thread
        // the block needs to be calculated first so iteration pointers can
        // be passed to be 

        thread_pool.emplace_back(
            compute_structures,
            args,
            args,
            args
        )
    }
    for (auto% t : thread_pool){
        t.join();
    }

}

void Gene::compute_blocked_structures(Model $model, structure_start, structure_stop, seq_start,){

    // the sequence is not directly connected to the structures which is
    // kind of annoying

    start = sequence.start() + *structure.position.start_pos
    end = sequence.end() + length of sequence which different for circular DNA

    model.compute_structure(sequence,start,stop,first_structure);

    model.compute_structure(*structure_start)




}

void Gene::create_structures(Model &model){

    std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
    windower.reset_window();
    int block_counter = 0

    windower.next_window_from_all_windows(start,stop);
    Structure first_structure;
    first_structure.position.chromosome = position.chromosome;
    first_structure.position.strand = position.strand;
    first_structure.position.start_pos = position.start_pos + windower.get_current_start_offset();
    first_structure.position.end_pos = position.start_pos + windower.get_current_stop_offset();

    rloop_structures.push_back(first_structure);

    // add the first structure to get going

    while (windower.has_next_window()){

        block_counter += windower.next_window_from_all_windows(start,stop);
        // will only be one when entering a new block

        Structure temp;
        temp.position.chromosome = position.chromosome;
        temp.position.strand = position.strand;
        temp.position.start_pos = position.start_pos + windower.get_current_start_offset();
        temp.position.end_pos = position.start_pos + windower.get_current_stop_offset();
        temp.block_id = block_counter;  // need to add this attribute to structure

        rloop_structures.push_back(temp);

        
       // cout << temp.position.start_pos << " " << temp.position.end_pos << endl;
        // do not compute the structure yet only create the structures
        // if (flag == 1){  // Structure is being built from new start site
        //     // is a parent structure so do things normally
        //     model.compute_structure(sequence,start,stop,temp);
        // }
        // else{
        //     // compute free energy and boltz using the previous structure
        //     model.compute_structure(sequence,start,stop,rloop_structures.back(), temp);
        // }      
}


void Gene::compute_structures(Model &model){

    // if positive strand is revesed at this point

    vector<char> temp_circular_sequence;
    if (sequence.size() == 0){
        //throw exception
    }
    //initializing the iterators ensures that the intial comparison in next_window_from_all_windows is not problematic
    std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
    // would need to do this for each thread
    windower.reset_window();
    // Get things rolling by doing first calculation outside of the loop    
    windower.next_window_from_all_windows(start,stop);
    Structure first_structure;
    first_structure.position.chromosome = position.chromosome;
    first_structure.position.strand = position.strand;
    first_structure.position.start_pos = position.start_pos + windower.get_current_start_offset();
    first_structure.position.end_pos = position.start_pos + windower.get_current_stop_offset();


    //cout << *start << endl;
    //cout << *stop << endl;

    model.compute_structure(sequence,start,stop,first_structure);

    rloop_structures.push_back(first_structure);

    int i = 0;
    while (windower.has_next_window()){

        int flag = windower.next_window_from_all_windows(start,stop);

        Structure temp;
        temp.position.chromosome = position.chromosome;
        temp.position.strand = position.strand;
        temp.position.start_pos = position.start_pos + windower.get_current_start_offset();
        temp.position.end_pos = position.start_pos + windower.get_current_stop_offset();
        
       // cout << temp.position.start_pos << " " << temp.position.end_pos << endl;
        
        if (flag == 1){  // Structure is being built from new start site
            // is a parent structure so do things normally
            model.compute_structure(sequence,start,stop,temp);
        }
        else{
            // compute free energy and boltz using the previous structure
            model.compute_structure(sequence,start,stop,rloop_structures.back(), temp);
        }
        // cout << "Got one here by guy " << i << endl;
        // cout << temp.free_energy << endl;
        // cout << temp.position.start_pos << endl;
        // cout << temp.position.end_pos << endl;
        // i ++;


        // here need to do the actual calculations by refering to the last
        // structures in some way
        // how do we want to do that by the way is the though thing

        // need to compute the first structure in this section of structs
        // and then can use that one to calculate the rest

        /*
        if current_boltz and current_free energy are 0:
            Means we started a new block (need signal for that)
            Create the first structure like normal
            Compute all other structures from the first one
        */
       rloop_structures.push_back(temp);

    }

    ground_state_energy = model.ground_state_energy();  // O(1)
}

void Gene::compute_structures_circular(Model &model){
    //check for circular sequence conditions
    if (sequence.size() == 0){
        //throw exception
		
    }
    //initializing the iterators ensures that the intial comparison in next_window_from_all_windows is not problematic
    std::vector<char>::iterator start = sequence.begin(),stop=sequence.begin()+1;
    windower.reset_window();
    while (windower.has_next_window_circular()){
        windower.next_window_from_all_windows_circular(start,stop);
        Structure temp;
        //set the Loci of the structure using the gene's Loci
        temp.position.chromosome = position.chromosome;
        temp.position.strand = position.strand;
        temp.position.start_pos = position.start_pos + windower.get_current_start_offset();
        temp.position.end_pos = position.start_pos + windower.get_current_stop_offset();
        //pass the structure and window boundaries to the model
        model.compute_structure(sequence,start,stop,temp);
        //push the now computed structure onto these_structures
        rloop_structures.push_back(temp);
    }
    //cout << rloop_structures.size() << endl;
}

void Gene::compute_residuals(Model &model){
    //verify that the structures have been computed
    //iterate through all the structures
    for (int i=0; i < rloop_structures.size(); i++){ //not iterating through all the structures???
        model.compute_residuals(rloop_structures[i]);
    }
}

void Gene::clear_structures(){
    rloop_structures.clear();
}

void Gene::complement_sequence(){
    for (int i = 0; i < sequence.size(); i++){
        if (sequence[i] == 'A')
            sequence[i] = 'T';
        else if (sequence[i] == 'T')
            sequence[i] = 'A';
        else if (sequence[i] == 'C')
            sequence[i] = 'G';
        else //sequence_data[i] == 'G'
            sequence[i] = 'C';
    }
}

void Gene::invert_sequence(){
    char temp;
    for (int i = 0; i < (sequence.size()/2); i++){
        temp = sequence[i];
        sequence[i] = sequence[sequence.size() - 1 - i];
        sequence[sequence.size() - 1 - i] = temp;
    }
}

int Gene::get_length(){
    return sequence.size();
}

void Gene::clear_sequence(){
    //delete sequence data
    sequence.clear();
}

void Gene::dump_structures(string outfilename){
    ofstream dumpfile(outfilename + gene_name + "_dump.txt",ios::out);
    std::stringstream ss;
    ss << "start_position stop_position energy probability\n";
    ss << "0 0 " << ground_state_energy << ' ' << 0 << endl; //add the ground state probabilit
    if (position.strand == "+") {
        for (int i = 0; i < rloop_structures.size(); i++) {
            ss << rloop_structures[i].position.start_pos << ' ' << rloop_structures[i].position.end_pos << ' ' <<
               rloop_structures[i].free_energy << ' ' << rloop_structures[i].probability << endl;
        }
    } else if (position.strand == "-") {
        for (int i = 0; i < rloop_structures.size(); i++) {
            ss << (this->getPosition().end_pos-rloop_structures[i].position.end_pos+this->getPosition().start_pos) << ' ' << (this->getPosition().end_pos-rloop_structures[i].position.start_pos+this->getPosition().start_pos) << ' ' <<
               rloop_structures[i].free_energy << ' ' << rloop_structures[i].probability << endl;
        }
    }
    else {
        cout << "Dump error. Strand unspecified.";
        exit(1); //replace with exception
    }
    dumpfile << ss.rdbuf();
    dumpfile.close();
}

