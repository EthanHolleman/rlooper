//
// Created by Robert Stolz on 6/28/17.
//

#include "windower.h"

Windower::Windower(): min_window_size(2), is_circular(false) {};

Windower::Windower(std::vector<char> &target_sequence){
    current_sequence = &target_sequence;
    current_start = target_sequence.begin();
    current_stop = current_start + min_window_size-1;
    is_circular = false;
    sliding_window_size = 0;
}

int Windower::get_min_window_size(){
    return min_window_size;
}

void Windower::set_min_window_size(int size){
    if (size < 2){
        throw WindowerException();
    }
    min_window_size = size;
}

void Windower::set_circular(bool Is_circular){
    is_circular = Is_circular;
}

void Windower::set_sequence(std::vector<char>& target_sequence){ //not working?
    current_sequence = &target_sequence;
    current_start = target_sequence.begin();
    current_stop = current_start + min_window_size - 2;
    sequence_size = target_sequence.size();
    current_window_size = min_window_size - 2;
}

bool Windower::has_next_window(){
    //if the last window has been reached
    if (is_circular){
        return has_next_window_circular();
    }
    if (current_start == current_stop-min_window_size+1 && current_stop == current_sequence->end()-1){
    
        /*
        NOTE: We do not have a next window when ...

        The window start is equal to the (window end - window size)

        -----|--------------------|------  The window

        min window size is in fact the min window size so windows can be larger
        which in term seems to effect iteration and window creation somewhere
        else (RUNTIME!!!)

        Is there some op for memoization / caching to some degree?
        Is there some connection between

        Window A: ABCD | EFGHIJKL | MNOP

        and 

        Window B: ABCDE | FGHIJKL | MNOP

        Or in otherwords does window FGHIJKL + E == EFGHIJKL

        If so then we can dynamic program the shit out of this but I'm not
        actually sure that is how things are working out

        *EH


        */



        return false;
    }
    else
        return true;
}

bool Windower::has_next_window_circular(){
    //if the last window has been reached
    if (current_start == current_sequence->end()-1 && current_window_size == sequence_size){
        return false;
    }
    else
        return true;
}

//change function type to return int *EH
int Windower::next_window_from_all_windows(std::vector<char>::iterator& start, std::vector<char>::iterator& stop){
    //if (!has_next_window()){ //safety check to make sure the next window exists
    //    throw WindowerException(); //throw exception?
    //}
    int flag = 0;
    if (is_circular){
        return next_window_from_all_windows_circular(start,stop);
    }
    if (stop < current_sequence->end()-1){  // Move the end further along but keep same start *EH
        ++current_stop;
        // flag that we have not yet moved to new block
    }
    else{ //if (start < current_sequence->end()-min_window_size){
        // We hit the current stop so we need to increment the start by one
        // And move the next stop by 1 min_window_size unit along the sequence
        ++current_start;
        current_stop = current_start + min_window_size-1;
        flag = 1;
    }
    start = current_start;
    stop = current_stop;
    return flag;
    //print_current_window();

    /*
    Think the above looks something like this

    min_window_size = 3
    start = 0 
    stop = 2 (min_window_size - 1)

    1. | A B C | D E

    Increment current_stop (stop)

    stop = 3

    2. | A B C D | E

    Increment current_stop (stop)

    stop = 4

    3. | A B C D E |

    Increment current start

    start = 1

    stop = start (1) + min_window_size (3) -1 = 3

    4. A | B C D | E

    ...

    *EH 

    */
}

int Windower::next_window_from_all_windows_circular(std::vector<char>::iterator& start, std::vector<char>::iterator& stop){
    //if (!has_next_window()){ //safety check to make sure the next window exists
    //    throw WindowerException(); //throw exception?
    //}
    int flag = 0;
    if (current_window_size < sequence_size){ //sequence length has not been reached
        if (current_stop == current_sequence->end()-1){ //boundary condition
            current_stop = current_sequence->begin();
        }
        else{
            ++current_stop;
        }
        current_window_size++;
    }
    else{ //implies sequence length has been reached, both start and stop need to be reassigned.
        ++current_start;
        current_stop = current_start + min_window_size-1;
        if (current_start == current_sequence->end()-1){ //if current start is the last position in the sequence
            current_stop = current_sequence->begin() + min_window_size-2;
        }
        current_window_size = min_window_size;
    }
    start = current_start;
    stop = current_stop;
    flag = 1;
    //print_current_window();
    return flag;
}

long int Windower::get_current_start_offset(){
    return current_start - current_sequence->begin();
}

long int Windower::get_current_stop_offset(){
    return current_stop - current_sequence->begin();
}

void Windower::reset_window(){
    current_start = current_sequence->begin();
    current_stop = current_sequence->begin()+min_window_size - 2;
    current_window_size = min_window_size-1;
}

void Windower::print_current_window(){
    std::vector<char>::iterator it;
    it = current_start;
    while (it != current_stop){
        //circular sequence boundary awareness condition here
        if (it == current_sequence->end()){ //boundary condition
            it = current_sequence->begin();
        }
        else {
            std::cout << *it << '\n';
            ++it;
        }
    }
    std::cout << *it << '\n';
    std::cout << "=====================\n";
}