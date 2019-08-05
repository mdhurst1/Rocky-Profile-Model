// SL.hpp

#ifndef SL_HPP
#define SL_HPP

class SL
{
    private:
        //internal variables
        double SeaLevel;

        // private initialise function
        void Initialise();
        void Initialise(double InputSL);

    public:
        // these are called constructors since it is a call to itself
        // we can have empty constructors or constructors with arguments
        SL() {Initialise();}
        SL(double InputSL) {Initialise(InputSL);}

        // function to get a value from inside the object
        double get_SeaLevel() {return SeaLevel;}
};

#endif