#ifndef Dco_hpp_
#define Dco_hpp_

#include "AlpsConfig.h"
#include "BcpsConfig.h"
#include "DcoConfig.hpp"

//! \page handle HomePage

/*! \mainpage

  # Bcps ideas for future

  We keep relaxed cols/rows (the whole object or integrality of the object) in
  BcpsModel in list BcpsObject ** relaxed_. This way Bcps can check the
  feasibility of its relaxed objects through virtual functions that will be
  impelemnted in Disco level.

  Can subproblems have different relaxed objects? Yes they can. But make sure
  all relaxed objects are in the relaxed_. If subproblems have different
  relaxed objects we will get a more depocposition like algorithm. Moreover if
  not all integrality is relaxed we might need a discrete solver (a DcoModel
  instance?) to solve the problem.

  Subproblems might have different solvers? A subprobllem might be an LP, SOCO,
  MILP or even MISOCO. How will DcoTreeNode::bound know about the solver to be
  used?

  When subproblem is solved with the solver we go back and check the
  feasibility of the objects that are relaxed in the subproblem (objects in
  relaxed_). This is done through infeasible() virtual function defined in
  BcpsObject.
  After this check we have a list of objects that are infeasible. At this
  point we need to decide what to do with them. Options are (1) generate
  cuts (using generateConstraints()), (2) lift the subproblem (using
  generateVariables()), (3) branch.
  (1) We can generate cuts when infeasible objects are continuous or
  integer. Generate a separating hyperplane that cuts the subproblem solution
  from feasible region.
  (2) Branching when the object is continuous. This is similar to branch
  on variables. Branching procedure should create new subproblems where
  infeasible object has new upper and lower bounds.

  feasibility checking should be implemented in Bcps level. Itertate over
  cols/rows and check their feasiblity. Store infeasible cols (BcpsVariable)
  and rows (BcpsConstraint). This function checks the feasibility of the
  solution stored in the solver at the time of the call.

  BcpsModel::feasibleSolution() should be able to check cols or rows only.
  In a typical branch and bound we need to check feasibility of cols only.
  In DisCO we may want to check both or check cols only.

  # DisCO ideas

  # Style Guide

  <ul>

    <li> Use two space indentation. Ted might object this, but this is my
    (aykut) favorite anyway. He will use commit hooks to get code in his
    favorite form.

    <li> Avoid too many nested code blocks and many indentation levels. Reading
    code with too many nested blocks is not fun. 3 levels of for/if blocks as
    follows is fine

    ```
    void ClassName::func() {
      for (;;) {
        if (boolean1) {
          for (;;) {
          }
        }
      }
    }
    ```

    If you have 4 or more, most probably there is a better way of doing what you
    are doing (create a new function, design change, etc.).

    <li> Curly braces that mark start of code block is used after the block
    specifier as in the example from the previous bullet. This is true for
    functions and class definitions too.

    <li> One line between function definitions. Two lines between class
    declerations if the header file has more than one.

    <li> One class decleration per header file. It is OK to add two classes in
    the same header if (1) they are closely related and (2) the second one is
    not a major class but just a helper for the major class (class that has the
    same name as the hearder file).

    <li> Class declared in a header file should have the same name as the
    header file. Source file that contains definition of the class should have
    the same name.

    <li> No whitespace at line endings. Please do not insert whitespaces at
    line endings. Such trailing white spaces should be removed. Trailing white
    spaces are redundant characters and are not wellcomed in DisCO source
    code. I (aykut) use following snippet in my .emacs, which removes the
    trailing white spaces at each save. I strongly advise you to do so.

    ```
    ;; remove trailing whitespace at each save
    (add-hook 'before-save-hook 'delete-trailing-whitespace)
    ```

    <li> If a function does not fit in the same line of its decleration in the
    header file then it is not supposed to be there. Move it to .cpp file.

    <li> All functions defined in .hpp are already inline since bla.bla version
    of gcc. This makes inline keywords in the header files redundant. Do not
    use inline for these functions, I (aykut) can not bear redundant stuff
    floating around.

    <li> Define pointers as `Type * name`. Asterisk character is neither next
    to `Type` nor variable `name`. `*` is separate since it is neither part of
    the `Type` nor variable `name`. This is the way recommended by Bjarne
    Stroustrup and it makes sense.

    <li> We use `const` specifiers as suggested by Bjarne Stroustrup in his
    book.  Check DisCO source code to see how const objects are declared, check
    Bjarne Stroustrup's "The C++ Language" book to see why.

    <li> Functions in public interface use
    <a>href="https://en.wikipedia.org/wiki/Camel_case">Camel Case</a> starting
    with a lower case. Same is true for class members.

    <li> From time to time private function might use underscores.

    <li> Right now we do not have a rule for local variables in
    functions. Camel Case and underscores are used together. This might
    change. Not sure about this now.

    <li> Header files should have "hpp" as extension. Source files should have
    "cpp" as extension.

    <li> Each header file should have include safeguard macros. Include
    safeguard have a specific format. For header file Dco.hpp include safeguard
    is as follows.

    ```
    #ifndef Dco_hpp_
    #define Dco_hpp_
    Stuff that is needed in header file.
    #endif
    ```

    <li> In source files limit rows to 79 characters. There might be a few
    lines ignoring this for now, this should be fixed. Following might be
    usefull in case you are using emacs (M-q wraps lines).

    ```
    (setq-default fill-column 79)
    ```

    <li> Adding comments for a specific part of code. To add comments for a
    specific part of code, first put the relevant code snippet into a block and
    write comments for the block. Comments should be in doxygen format. See
    following example

    ```
    start doxygen comments (forward slash, asterisk, exclamation)
     Comments regarding the following code block.
    end doxygen comments (asterisk, forward slash)
    {
      code that is being commented
    }
    ```

    <li> Go local as much as possible. Create variables in the namespace level
    that is directly relevant not in the upper levels.

    <li> Keep documentation up to date.

  </ul>

*/

#define DISCO_INFINITY 1e30

enum DcoConstraintType {
  DcoConstraintTypeLinear = 0,
  DcoConstraintTypeConic
};

enum DcoLorentzConeType {
  DcoLorentzCone = 0,
  DcoRotatedLorentzCone
};

enum DcoSubproblemStatus{
  DcoSubproblemStatusOptimal,
  DcoSubproblemStatusAbandoned,
  DcoSubproblemStatusPrimalInfeasible,
  DcoSubproblemStatusDualInfeasible,
  DcoSubproblemStatusPrimalObjLim,
  DcoSubproblemStatusDualObjLim,
  DcoSubproblemStatusIterLim,
  DcoSubproblemStatusUnknown
};

enum DcoReturnStatus {
  DcoReturnStatusOk = 0,
  DcoReturnStatusErrLp,
  DcoReturnStatusInfeasible,
  DcoReturnStatusUnbounded,
  DcoReturnStatusOverObjLim,
  DcoReturnStatusFeasible,
  DcoReturnStatusBranch,
  DcoReturnStatusUnknown
};


enum DcoCutStrategy {
  DcoCutStrategyNotSet = -1,
  DcoCutStrategyNone = 0,
  DcoCutStrategyRoot,
  DcoCutStrategyAuto,
  DcoCutStrategyPeriodic
};

enum DcoConicCutStrategy {
  DcoConicCutStrategyNotSet = -1,
  DcoConicCutStrategyNone = 0,
  DcoConicCutStrategyRoot,
  DcoConicCutStrategyAuto,
  DcoConicCutStrategyPeriodic
};

enum DcoHeurStrategy {
  DcoHeurStrategyNotSet = -1,
  DcoHeurStrategyNone = 0,
  DcoHeurStrategyRoot,
  DcoHeurStrategyAuto,
  DcoHeurStrategyPeriodic,
  DcoHeurStrategyBeforeRoot // Before solving first relaxation
};

enum DcoHeurType {
  DcoHeurTypeNotSet = -1,
  DcoHeurTypeRounding
};

enum DcoHotStartStrategy{
  DcoHotStartBranchIncorrect,
  DcoHotStartBranchCorrect
};

enum DcoBranchingStrategy{
  DcoBranchingStrategyMaxInfeasibility,
  DcoBranchingStrategyPseudoCost,
  DcoBranchingStrategyReliability,
  DcoBranchingStrategyStrong,
  DcoBranchingStrategyBilevel
};

enum DcoSolutionType {
  DcoSolutionTypeBounding,
  DcoSolutionTypeBranching,
  DcoSolutionTypeDiving,
  DcoSolutionTypeHeuristic,
  DcoSolutionTypeStrong
};

/** Branching object type. */
enum DcoBranchingObjectType {
  DcoBranchingObjectTypeNone = 0,
  DcoBranchingObjectTypeInt,
  DcoBranchingObjectTypeSos,
  DcoBranchingObjectTypeBilevel
};

/** Node branch direction, is it a left node or right */
enum DcoNodeBranchDir {
  DcoNodeBranchDirectionDown = 0,
  DcoNodeBranchDirectionUp
};

/** Integral type */
enum DcoIntegralityType {
  DcoIntegralityTypeCont = 0,
  DcoIntegralityTypeInt
};

#endif
