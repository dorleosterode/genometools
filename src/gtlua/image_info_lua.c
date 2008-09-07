/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef WITHOUT_CAIRO

#include <float.h>
#include <limits.h>
#include "lauxlib.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/recmap.h"
#include "core/error.h"
#include "extended/luahelper.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/image_info_lua.h"

static int imageinfo_lua_new(lua_State *L)
{
  GT_ImageInfo **ii;
  ii = lua_newuserdata(L, sizeof (GT_ImageInfo*));
  assert(ii);
  *ii = gt_image_info_new();
  luaL_getmetatable(L, IMAGEINFO_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int imageinfo_lua_get_height(lua_State *L)
{
  GT_ImageInfo **ii;
  unsigned long height;
  GT_Error *err = error_new();
  ii = check_imageinfo(L, 1);
  assert(ii);
  height = gt_image_info_get_height(*ii);
  if (height > DBL_MAX)
  {
    error_set(err, "image height exceeds %f!", DBL_MAX);
    return lua_gt_error(L, err);
  }
  else
    lua_pushnumber(L, (double) height);
  error_delete(err);
  return 0;
}

static int imageinfo_lua_num_of_recmaps(lua_State *L)
{
  GT_ImageInfo **ii;
  unsigned long nof_rm;
  GT_Error *err = error_new();
  ii = check_imageinfo(L, 1);
  assert(ii);
  nof_rm = gt_image_info_num_of_recmaps(*ii);
  if (nof_rm > DBL_MAX)
  {
    error_set(err, "number of recmaps exceeds %f!", DBL_MAX);
    return lua_gt_error(L, err);
  }
  else
    lua_pushnumber(L, (double) nof_rm);
  error_delete(err);
  return 0;
}

static void push_recmap_as_table(lua_State *L, const GT_RecMap *rm)
{
  assert(rm);
  lua_newtable(L);
  lua_pushstring(L, "nw_x");
  lua_pushnumber(L, rm->nw_x);
  lua_rawset(L, -3);
  lua_pushstring(L, "nw_y");
  lua_pushnumber(L, rm->nw_y);
  lua_rawset(L, -3);
  lua_pushstring(L, "se_x");
  lua_pushnumber(L, rm->se_x);
  lua_rawset(L, -3);
  lua_pushstring(L, "se_y");
  lua_pushnumber(L, rm->se_y);
  lua_rawset(L, -3);
  lua_pushstring(L, "feature_ref");
  genome_node_lua_push(L, genome_node_rec_ref((GenomeNode*) rm->gf));
  lua_rawset(L, -3);
}

static int imageinfo_lua_recmaps_as_table(lua_State *L)
{
  GT_ImageInfo **ii;
  unsigned long num, i;
  ii = check_imageinfo(L, 1);
  assert(ii);
  num = gt_image_info_num_of_recmaps(*ii);
  if (num>0)
  {
    lua_newtable(L);
    for (i=0;i<num;i++)
    {
      lua_pushnumber(L, i+1);
      push_recmap_as_table(L, gt_image_info_get_recmap(*ii, i));
      lua_rawset(L, -3);
    }
  } else lua_pushnil(L);
  return 1;
}

static int imageinfo_lua_delete(lua_State *L)
{
  GT_ImageInfo **ii;
  ii = check_imageinfo(L, 1);
  gt_image_info_delete(*ii);
  return 0;
}

static const struct luaL_Reg imageinfo_lib_f [] = {
  { "imageinfo_new", imageinfo_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg imageinfo_lib_m [] = {
  { "get_height", imageinfo_lua_get_height },
  { "num_of_recmaps", imageinfo_lua_num_of_recmaps },
  { "get_recmaps", imageinfo_lua_recmaps_as_table },
  { NULL, NULL }
};

int luaopen_imageinfo(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, IMAGEINFO_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, imageinfo_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, imageinfo_lib_m);
  luaL_register(L, "gt", imageinfo_lib_f);
  return 1;
}

#endif
