# echo "cleanup TestAlgo TestAlgo-00-00-01 in /afs/ihep.ac.cn/users/p/privalov/batch/7.0.3"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTestAlgotempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTestAlgotempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  $* >${cmtTestAlgotempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  $* >${cmtTestAlgotempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtTestAlgotempfile}
  unset cmtTestAlgotempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtTestAlgotempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtTestAlgotempfile}
unset cmtTestAlgotempfile
return $cmtcleanupstatus

