// RP.hpp

#ifndef RP_HPP
#define RP_HPP

class RP
{
    private:
        //internal variables
        double Gradient;

        // private initialise functions
        void Initialise();
        void Initialise(double InputGradient);

    public:
        // these are called constructors since it is a call to itself
        // we can have empty constructors or constructors with arguments
        RP() {Initialise();}
        RP(double InputGradient) {Initialise(InputGradient);}

        // function to get a value from inside the object
        double get_Gradient() {return Gradient;}
};

#endif