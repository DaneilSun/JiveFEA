
program = solid

subdirs = JJ_Material JJ_Models JJ_Modules JJ_Util JJ_Fixes

include $(JEMDIR)/makefiles/packages/base.mk

include $(JIVEDIR)/makefiles/packages/algebra.mk
include $(JIVEDIR)/makefiles/packages/app.mk
include $(JIVEDIR)/makefiles/packages/fem.mk
include $(JIVEDIR)/makefiles/packages/geom.mk
include $(JIVEDIR)/makefiles/packages/gl.mk
include $(JIVEDIR)/makefiles/packages/implict.mk
include $(JIVEDIR)/makefiles/packages/model.mk
include $(JIVEDIR)/makefiles/packages/solver.mk
include $(JIVEDIR)/makefiles/packages/util.mk

include $(JIVEDIR)/makefiles/prog.mk

MY_INCDIRS = . $(subdirs)
