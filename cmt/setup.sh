# echo "setup TestAlgo TestAlgo-00-00-01 in /afs/ihep.ac.cn/users/p/privalov/batch/7.0.3"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTestAlgotempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTestAlgotempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  -no_cleanup $* >${cmtTestAlgotempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  -no_cleanup $* >${cmtTestAlgotempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtTestAlgotempfile}
  unset cmtTestAlgotempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtTestAlgotempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
export LD_LIBRARY_PATH=$TESTALGOROOT/x86_64-slc6-gcc46-opt/:$LD_LIBRARY_PATH
/bin/rm -f ${cmtTestAlgotempfile}
unset cmtTestAlgotempfile
return $cmtsetupstatus

