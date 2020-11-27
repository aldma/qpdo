#ifndef NEWTON_H
    #define NEWTON_H

    # ifdef __cplusplus
        extern "C" {
    # endif

    #include "types.h"

    void newton_direction(QPDOWorkspace *work, cholmod_common *c);

    void active_constraints(QPDOWorkspace *work);

    void enter_leave_constraints(QPDOWorkspace *work);

    # ifdef __cplusplus
        }
    # endif

#endif