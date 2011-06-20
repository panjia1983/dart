#ifndef OPTIMIZER_PROBLEM_H
#define OPTIMIZER_PROBLEM_H

#include <vector>

namespace optimizer {
    class Var;
    class ConstraintBox;
    class ObjectiveBox;
    
    class Problem {
    public:
        Problem();
        virtual ~Problem();
        virtual void update() = 0;

        void AddVariable(double value, double lower, double upper);
        void createBoxes();

    private:
        ConstraintBox* mConBox;
        ObjectiveBox* mObjBox;
        std::vector<Var *> mVariables;
        
    };
} // namespace optimizer

#endif // #ifndef OPTIMIZER_PROBLEM_H
